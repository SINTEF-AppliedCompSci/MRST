function [problem, state] = equationsBlackOilSolvent(state0, state, model, dt, drivingForces, varargin)
% Generate linearized problem for the black-oil solvent model equations
%
% SYNOPSIS:
%   [problem, state] = equationsBlackOilSolvent(state0, state, model, dt, drivingForces)
%
% DESCRIPTION:
%   This is the core function of the black-oil solvent model. This function
%   assembles the residual equations for the conservation of water, oil and
%   gas, as well as required well equations. By default, Jacobians are also
%   provided by the use of automatic differentiation.
%
%   The model assumes two extrema: In cases with only oil, reservoir gas
%   and water, we have traditional black-oil behavior. In cases with oil,
%   solvent gas and water, the oil mixes with the solvent gas according to
%   the Todd-Longstaff model. In intermediate regions, we interpolate
%   between the two extrema based on the solvet to total gas saturation
%   fraction, and the pressure.
%
% REQUIRED PARAMETERS:
%   state0    - Reservoir state at the previous timestep. Assumed to have
%               physically reasonable values.
%
%   state     - State at the current nonlinear iteration. The values do not
%               need to be physically reasonable.
%
%   model     - BlackOilSolventModel-derived class. Typically,
%               equationsBlackOil will be called from the class
%               getEquations member function.
%
%   dt        - Scalar timestep in seconds.
%
%   drivingForces - Struct with fields:
%                   * W for wells. Can be empty for no wells.
%                   * NOTE: The current implementation does not support
%                   BC's and sources.
%
% OPTIONAL PARAMETERS:
%   'Verbose'    -  Extra output if requested.
%
%   'reverseMode'- Boolean indicating if we are in reverse mode, i.e.
%                  solving the adjoint equations. Defaults to false.
%
%   'resOnly'    - Only assemble residual equations, do not assemble the
%                  Jacobians. Can save some assembly time if only the
%                  values are required.
%
%   'iterations' - Nonlinear iteration number. Special logic happens in the
%                  wells if it is the first iteration.
% RETURNS:
%   problem - LinearizedProblemAD class instance, containing the water, oil
%             and gas conservation equations, as well as well equations
%             specified by the FacilityModel class.
%
%   state   - Updated state. Primarily returned to handle changing well
%             controls from the well model.
%
% SEE ALSO:
%   BlackOilSolventModel,  ThreePhaseBlackOilModel


%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
                                                 
opt = struct('Verbose'    , mrstVerbose, ...
             'reverseMode', false      ,...
             'resOnly'    , false      ,...
             'iteration'  , -1         );

opt = merge_options(opt, varargin{:});

W     = drivingForces.W;
op    = model.operators;
fluid = model.fluid;

% Properties at current timestep
[p , sW , sG , sS , rs , rv , wellSol ] = model.getProps(state, ...
             'pressure', 'water', 'gas', 'solvent', 'rs', 'rv', 'wellSol');
% Properties at previous timestep
[p0, sW0, sG0, sS0, rs0, rv0, wellSol0] = model.getProps(state0, ...
             'pressure', 'water', 'gas', 'solvent', 'rs', 'rv', 'wellSol');

% Well variables
[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

% Initialize primary variables (gas saturation taken to be sG + sS)
st  = model.getCellStatusVO(state,  1-(sW + sG  + sS) , sW  + sS , sG );
st0 = model.getCellStatusVO(state0, 1-(sW0+ sG0 + sS0), sW0 + sS0, sG0);

if model.disgas || model.vapoil
    % X is either Rs, Rv or Sg, depending on each cell's saturation status
    x = st{1}.*rs + st{2}.*rv + st{3}.*sG;
    gvar = 'x';
else
    x = sG;
    gvar = 'sG';
end

if ~opt.resOnly,
    if ~opt.reverseMode,
        % define primary varible x and initialize
        [p, sW, x, sS, wellVars{:}] = initVariablesADI(p, sW, x, sS, wellVars{:});
    else
        wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
        x0 = st0{1}.*rs0 + st0{2}.*rv0 + st0{3}.*sG0;
        [p0, sW0, x0, sS0, wellVars0{:}] = initVariablesADI(p0, sW0, x0, sS0, wellVars0{:}); %#ok
        [sG0, rs0, rv0] = calculateHydrocarbonsFromStatusBO(model, st0, 1-(sW0 + sS0), x0, rs0, rv0, p0);
    end
end
    
if ~opt.reverseMode
    % Compute values from status flags. If we are in reverse mode, these
    % values have already converged in the forward simulation.
    [sG, rs, rv, rsSat, rvSat] = calculateHydrocarbonsFromStatusBO(model, st, 1 - (sW + sS), x, rs, rv, p);
end

% We will solve for pressure, and water/gas/solvent saturations (oil
% saturation follows via the definition of saturations), and well rates +
% bhp.
primaryVars = {'pressure', 'sW', gvar, 'sS', wellVarNames{:}};

sO  = 1 - (sW  + sG  + sS );
sO0 = 1 - (sW0 + sG0 + sS0);

% Calculate residual saturations
[sWcon, sOr , sGc ] = computeResidualSaturations(model, p , sW , sG , sS , state0);
if isfield(state0, 'sOr') && isfield(state0, 'sGcr')
    sOr0 = state0.sOr;
    sGc0 = state0.sGc;
else
    [~    , sOr0, sGc0] = computeResidualSaturations(model, p0, sW0, sG0, sS0, state0);
end

% Get multipliers
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);
T = op.T.*transMult;

 % Calculate effective relperms
[krW, krO, krG, krS] = computeRelPermSolvent(model, p, sW, sO, sG, sS, sWcon, sOr, sGc, mobMult);

 % Calulate effective viscosities and densities
[muW, muO, muG, muS, rhoW, rhoO , rhoG , rhoS , bW , bO , bG , bS , pW, pG] ...
    = computeViscositiesAndDensities(model, p , sW, sO , sG , sS , sOr , sGc , rs, rv, ~st{1}, ~st{2});
[~  , ~  , ~  , ~  , ~   , rhoO0, rhoG0, rhoS0, bW0, bO0, bG0, bS0, ~ , ~] ...
    = computeViscositiesAndDensities(model, p0, sW0, sO0, sG0, sS0, sOr0, sGc0, rs0, rv0, ~st0{1}, ~st0{2});

gdz = model.getGravityGradient();

[vW, vO, vG, vS, mobW, mobO, mobG, mobS, upcW, upcO, upcG, upcS] ...
    = getFluxAndPropsSolvent(fluid, p, krW, krO, krG, krS, muW, muO, muG, muS, rhoW, rhoO, rhoG, rhoS, sW, sG, sS, T, gdz, op);

if model.outputFluxes
    state = model.storeFluxes(state, vW, vO, vG, vS);
end
if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, bG, bS);
    state = model.storeDensity(state, rhoW, rhoO, rhoG, rhoS);
    state = model.storeMobilities(state, mobW, mobO, mobG, mobS);
    state = model.storeUpstreamIndices(state, upcW, upcO, upcG, upcS);
end

% EQUATIONS ---------------------------------------------------------------
% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.
bWvW = op.faceUpstr(upcW, bW).*vW;
bOvO = op.faceUpstr(upcO, bO).*vO;
bGvG = op.faceUpstr(upcG, bG).*vG;
bSvS = op.faceUpstr(upcS, bS).*vS;

% Conservation of mass for water
water = (op.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + op.Div(bWvW);

% Conservation of mass for oil
if model.vapoil
    % The model allows oil to vaporize into the gas phase. The conservation
    % equation for oil must then include the fraction present in the gas
    % phase.
    rvbGvG = op.faceUpstr(upcG, rv).*bGvG;
    % Final equation
    oil = (op.pv/dt).*( pvMult .*(bO.* sO  + rv .*bG .*sG )-    ...
                        pvMult0.*(bO0.*sO0 + rv0.*bG0.*sG0) ) + ...
           op.Div(bOvO + rvbGvG);
else
    oil = (op.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) + op.Div(bOvO);
end

% Conservation of mass for gas
if model.disgas
    % The gas transported in the oil phase.
    rsbOvO = op.faceUpstr(upcO, rs).*bOvO;
    
    gas = (op.pv/dt).*( pvMult.* (bG.* sG  + rs.* bO.* sO) - ...
        pvMult0.*(bG0.*sG0 + rs0.*bO0.*sO0 ) ) + ...
        op.Div(bGvG + rsbOvO);
else
    gas = (op.pv/dt).*( pvMult.*bG.*sG - pvMult0.*bG0.*sG0 ) + op.Div(bGvG);
end

% Conservation of mass for solvent
solvent = (op.pv/dt).*( pvMult.*bS.*sS - pvMult0.*bS0.*sS0 ) + op.Div(bSvS);

state.sr = [value(sOr), value(sGc)];
state.sOr = value(sOr);
state.sGc = value(sGc);


eqs   = {water, oil, gas, solvent};
names = {'water', 'oil', 'gas', 'solvent'};
types = {'cell', 'cell', 'cell', 'cell'};

rho = {rhoW, rhoO, rhoG, rhoS};
mob = {mobW, mobO, mobG, mobS};
sat = {sW, sO, sG, sS};
dissolved = model.getDissolutionMatrix(rs, rv);

% Finally, add in and setup well equations
if ~isempty(W)
    [eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, types, wellSol0, wellSol, wellVars, wellMap, p, mob, rho, dissolved, {}, dt, opt);
end

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

end
