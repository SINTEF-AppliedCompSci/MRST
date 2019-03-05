function [problem, state] = equationsThreePhaseBlackOilSurfactant(state0, state, model, ...
    dt, drivingForces, varargin)
%
% SYNOPSIS:
%   function [problem, state] = equationsThreePhaseBlackOilSurfactant(state0, state, model, dt, drivingForces, varargin)
%
% DESCRIPTION:
%   Assemble the linearized equations for a blackoil-surfactant system,
%   computing both the residuals and the Jacobians. Returns the result as
%   an instance of the class LinearizedProblem which can be solved using
%   instances of LinearSolverAD.
%
%   A description of the modeling equations can be found in the directory
%   ad-eor/docs.
%
%
% PARAMETERS:
%   state0        - State at previous times-step
%   state         - State at current time-step
%   model         - Model instance
%   dt            - time-step
%   drivingForces - Driving forces (boundary conditions, wells, ...)
%   varargin      - optional parameters
%
% RETURNS:
%   problem - Instance of LinearizedProblem
%   state   - Updated state variable (fluxes, mobilities and more can be
%             stored, the wellSol structure is also updated in case of control switching)
%
% EXAMPLE:
%
% SEE ALSO: LinearizedProblem, LinearSolverAD, OilWaterSurfactantModel
%
%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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

opt = struct('Verbose', mrstVerbose, ...
    'reverseMode', false, ...
    'resOnly', false, ...
    'iteration', -1 );
opt = merge_options(opt, varargin{:});

%% Shorter names for some commonly used parts of the model and forces.
W     = drivingForces.W;
fluid = model.fluid;
op    = model.operators;
G     = model.G;

% Properties at current timestep
[p, sW, sG, rs, rv, c, cmax, wellSol] = model.getProps(state, ...
    'pressure', 'water', 'gas', 'rs', 'rv', 'surfactant', 'surfactantmax', 'wellsol');

% Properties at previous timestep
[p0, sW0, sG0, rs0, rv0, c0, cmax0, wellSol0] = model.getProps(state0, ...
    'pressure', 'water', 'gas', 'rs', 'rv', 'surfactant', 'surfactantmax', 'wellsol');

[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);
% Typically the primary well variables are :
%  - the phase well rates (qWell)
%  - the bottom hole pressures
%  - the surfactant concentrations, at injection and production wells,
%    contained in the wellVars, wellVarNames, wellMap structures

%% Initialization of primary variables ----------------------------------
st  = model.getCellStatusVO(state,  1-sW-sG,   sW,  sG);
st0 = model.getCellStatusVO(state0, 1-sW0-sG0, sW0, sG0);

if model.disgas || model.vapoil
    % X is either Rs, Rv or Sg, depending on each cell's saturation status
    x = st{1}.*rs + st{2}.*rv + st{3}.*sG;
    gvar = 'x';
else
    x = sG;
    gvar = 'sG';
end

% Initialize independent variables.
if ~opt.resOnly
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode
        [p, sW, x, c, wellVars{:}] = initVariablesADI(p, sW, x, c, wellVars{:});
    else
        x0 = st0{1}.*rs0 + st0{2}.*rv0 + st0{3}.*sG0;
        % Set initial gradient to zero
        wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
        [p0, sW0, x0, c0, wellVars0{:}] = initVariablesADI(p0, sW0, x0, c0, wellVars0{:}); %#ok
        clear zw;
        [sG0, rs0, rv0] = calculateHydrocarbonsFromStatusBO(model, st0, 1-sW, x0, rs0, rv0, p0);
    end
end

if ~opt.reverseMode
    % Compute values from status flags. If we are in reverse mode, these
    % values have already converged in the forward simulation.
    [sG, rs, rv, rsSat, rvSat] = calculateHydrocarbonsFromStatusBO(model, st, 1-sW, x, rs, rv, p);
end

% We will solve for pressure, water and gas saturation (oil saturation follows via
% the definition of saturations), surfactant concentration and well rates +
% bhp.
primaryVars = {'pressure', 'sW', gvar, 'surfactant', wellVarNames{:}};



%% EQUATIONS ---------------------------------------------------------------
pBH = wellVars{wellMap.isBHP};
% Compute fluxes and other properties for oil and water.
[dp, mob, upc, b, rho, pvMult, b0, pvMult0, T] = ...
    computeFluxAndPropsThreePhaseBlackOilSurfactant(model, p0, p, sW, sG, c, pBH, W, rs, rv, st);

% divide to water/surfactant-oil-gas three parts

dpW  = dp{1} ; dpO  = dp{2};
mobW = mob{1}; mobO = mob{2};
rhoW = rho{1}; rhoO = rho{2};
upcW = upc{1}; upcO = upc{2};
bW   = b{1}  ; bO   = b{2};
bW0  = b0{1} ; bO0  = b0{2};
mobSft =  mobW.*c;

% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.
vO     = -op.faceUpstr(upcO, mobO).*T.*dpO;
vW     = -op.faceUpstr(upcW, mobW).*T.*dpW;
vG     = -op.faceUpstr(upcG, mobG).*T.*dpG;
vSft   = -op.faceUpstr(upcW, mobSft).*T.*dpW;
bOvO   =  op.faceUpstr(upcO, bO).*vO;
bWvW   =  op.faceUpstr(upcW, bW).*vW;
bGvG   =  op.faceUpstr(upcg, bG).*vG;
bWvSft =  op.faceUpstr(upcW, bW).*vSft;


%% Conservation of mass for water
water = (op.pv/dt).*(pvMult.*bW.*sW - pvMult0.*bW0.*sW0) + op.Div(bWvW);

%% Conservation of mass for oil
sO  = 1 - sW - sG;
sO0 = 1 - sW0 - sG0;
if model.vapoil
    % The model allows oil to vaporize into the gas phase. The conservation
    % equation for oil must then include the fraction present in the gas
    % phase.
    rvbGvG = op.faceUpstr(upcg, rv).*bGvG;
    % Final equation
    oil = (op.pv/dt).*( pvMult.* (bO.* sO  + rv.* bG.* sG) - ...
        pvMult0.*(bO0.*sO0 + rv0.*bG0.*sG0) ) + ...
        op.Div(bOvO + rvbGvG);
else
    oil = (op.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) + op.Div(bOvO);
end

%% Conservation of mass for gas
if model.disgas
    % The gas transported in the oil phase.
    rsbOvO = og.faceUpstr(upco, rs).*bOvO;

    gas = (og.pv/dt).*( pvMult.* (bG.* sG  + rs.* bO.* sO) - ...
        pvMult0.*(bG0.*sG0 + rs0.*bO0.*sO0 ) ) + ...
        og.Div(bGvG + rsbOvO);
else
    gas = (og.pv/dt).*( pvMult.*bG.*sG - pvMult0.*bG0.*sG0 ) + og.Div(bGvG);
end

%% Computation of adsoprtion term
poro = model.rock.poro;
ads  = effads(c, cmax, fluid);
ads0 = effads(c0, cmax0, fluid);
ads_term = fluid.rhoRSft.*((1-poro)./poro).*(ads - ads0);

%% Conservation of surfactant in water:
surfactant = (op.pv/dt).*((pvMult.*bW.*sW.*c - pvMult0.*bW0.*sW0.*c0) +  ads_term) + ...
    op.Div(bWvSft);

%%

if model.extraStateOutput
    sigma = fluid.ift(c);
end

eqs   = {water, oil, gas, surfactant};
names = {'water', 'oil', 'gas', 'surfactant'};
types = {'cell', 'cell', 'cell', 'cell'};

rho = {rhoW, rhoO, rhoG};
mob = {mobW, mobO, mobG};
sat = {sW, sO, sG};
dissolved = model.getDissolutionMatrix(rs, rv);

[eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
    {pW, p, pG}, sat, mob, rho, ...
    dissolved, {c}, ...
    drivingForces);
% Finally, add in and setup well equations
[eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, ...
    names, types, wellSol0, ...
    wellSol, ...
    wellVars, wellMap, ...
    p, mob, rho, dissolved, ...
    {c}, dt, opt);
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

end

% Effective adsorption, depending of desorption or not
function y = effads(c, cmax, fluid)
if fluid.adsInxSft == 2
    y = fluid.surfads(max(c, cmax));
else
    y = fluid.surfads(c);
end
end
