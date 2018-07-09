function [problem, state] = equationsBlackOil(state0, state, model, dt, drivingForces, varargin)
% Generate linearized problem for the black-oil equations
%
% SYNOPSIS:
%   [problem, state] = equationsBlackOil(state0, state, model, dt, drivingForces)
%
% DESCRIPTION:
%   This is the core function of the black-oil solver. This function
%   assembles the residual equations for the conservation of water, oil and
%   gas, as well as required well equations. By default, Jacobians are also
%   provided by the use of automatic differentiation.
%
%   Oil can be vaporized into the gas phase if the model has the vapoil
%   property enabled. Analogously, if the disgas property is enabled, gas
%   is allowed to dissolve into the oil phase. Note that the fluid
%   functions change depending on vapoil/disgas being active and may have
%   to be updated when the property is changed in order to run a successful
%   simulation.
%
% REQUIRED PARAMETERS:
%   state0    - Reservoir state at the previous timestep. Assumed to have
%               physically reasonable values.
%
%   state     - State at the current nonlinear iteration. The values do not
%               need to be physically reasonable.
%
%   model     - ThreePhaseBlackOilModel-derived class. Typically,
%               equationsBlackOil will be called from the class
%               getEquations member function.
%
%   dt        - Scalar timestep in seconds.
%
%   drivingForces - Struct with fields:
%                   * W for wells. Can be empty for no wells.
%                   * bc for boundary conditions. Can be empty for no bc.
%                   * src for source terms. Can be empty for no sources.
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
%   equationsOilWater, ThreePhaseBlackOilModel

%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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
opt = struct('Verbose',     mrstVerbose,...
            'reverseMode', false,...
            'resOnly',     false,...
            'iteration',   -1);

opt = merge_options(opt, varargin{:});

% Shorter names for some commonly used parts of the model and forces.
s = model.operators;
f = model.fluid;

% Properties at current timestep
[p, sW, sG, rs, rv, wellSol] = model.getProps(state, ...
    'pressure', 'water', 'gas', 'rs', 'rv', 'wellSol');
% Properties at previous timestep
[p0, sW0, sG0, rs0, rv0, wellSol0] = model.getProps(state0, ...
    'pressure', 'water', 'gas', 'rs', 'rv', 'wellSol');

[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

if ~model.water
    [sW, sW0] = deal(0);
end

%Initialization of primary variables ----------------------------------
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

if ~opt.resOnly
    if ~opt.reverseMode
        % define primary varible x and initialize
        if model.water
            [p, sW, x, wellVars{:}] = model.AutoDiffBackend.initVariablesAD(p, sW, x, wellVars{:});
        else
            [p, x, wellVars{:}] = model.AutoDiffBackend.initVariablesAD(p, x, wellVars{:});
        end
    else
        wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
        x0 = st0{1}.*rs0 + st0{2}.*rv0 + st0{3}.*sG0;
        if model.water
            [p0, sW0, x0, wellVars0{:}] = model.AutoDiffBackend.initVariablesAD(p0, sW0, x0, wellVars0{:}); %#ok
        else
            [p0, x0, wellVars0{:}] = model.AutoDiffBackend.initVariablesAD(p0, x0, wellVars0{:}); %#ok
        end
        [sG0, rs0, rv0] = calculateHydrocarbonsFromStatusBO(model, st0, 1-sW0, x0, rs0, rv0, p0);
    end
end
% We will solve for pressure, water and gas saturation (oil saturation
% follows via the definition of saturations) and well rates + bhp.
if model.water
    primaryVars = {'pressure', 'sW', gvar, wellVarNames{:}};
else
    primaryVars = {'pressure', gvar, wellVarNames{:}};
end
    
if ~opt.reverseMode
    % Compute values from status flags. If we are in reverse mode, these
    % values have already converged in the forward simulation.
    [sG, rs, rv, rsSat, rvSat] = calculateHydrocarbonsFromStatusBO(model, st, 1-sW, x, rs, rv, p);
end
% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

% Evaluate relative permeability
sO  = 1 - sW  - sG;
sO0 = 1 - sW0 - sG0;
if model.water
    sat = {sW, sO, sG};
    [krW, krO, krG] = model.evaluateRelPerm(sat);
    % Mobility multiplier for water
    krW = mobMult.*krW;
else
    sat = {sO, sG};
    [krO, krG] = model.evaluateRelPerm(sat);
end

% Modifiy relperm by mobility multiplier (if any)
krO = mobMult.*krO; krG = mobMult.*krG;

% Compute transmissibility
T = s.T.*transMult;

% Gravity gradient per face
gdz = model.getGravityGradient();

% Evaluate water properties
if model.water
    [vW, bW, mobW, rhoW, pW, upcw] = getFluxAndPropsWater_BO(model, p, sW, krW, T, gdz);
    bW0 = f.bW(p0);
else
    [vW, bW, mobW, rhoW, pW, upcw] = deal([]);
end

% Evaluate oil properties
[vO, bO, mobO, rhoO, p, upco] = getFluxAndPropsOil_BO(model, p, sO, krO, T, gdz, rs, ~st{1});
bO0 = getbO_BO(model, p0, rs0, ~st0{1});

% Evaluate gas properties
bG0 = getbG_BO(model, p0, rv0, ~st0{2});
[vG, bG, mobG, rhoG, pG, upcg] = getFluxAndPropsGas_BO(model, p, sG, krG, T, gdz, rv, ~st{2});

% Store fluxes / properties for debugging / plotting, if requested.
if model.outputFluxes
    state = model.storeFluxes(state, vW, vO, vG);
end
if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, bG);
    state = model.storeDensity(state, rhoW, rhoO, rhoG);
    state = model.storeMobilities(state, mobW, mobO, mobG);
    state = model.storeUpstreamIndices(state, upcw, upco, upcg);
end

% EQUATIONS -----------------------------------------------------------

% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.
bOvO = s.faceUpstr(upco, bO).*vO;
bGvG = s.faceUpstr(upcg, bG).*vG;

% The first equation is the conservation of the water phase. This equation is
% straightforward, as water is assumed to remain in the aqua phase in the
% black oil model.
if model.water
    bWvW = s.faceUpstr(upcw, bW).*vW;
    water = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 );
    divWater = s.Div(bWvW);
end

% Second equation: mass conservation equation for the oil phase at surface
% conditions. This is any liquid oil at reservoir conditions, as well as
% any oil dissolved into the gas phase (if the model has vapoil enabled).
if model.vapoil
    % The model allows oil to vaporize into the gas phase. The conservation
    % equation for oil must then include the fraction present in the gas
    % phase.
    rvbGvG = s.faceUpstr(upcg, rv).*bGvG;
    % Final equation
    oil = (s.pv/dt).*( pvMult .* (bO.* sO  + rv.* bG.* sG) - ...
                       pvMult0.*(bO0.*sO0 + rv0.*bG0.*sG0));
    divOil = s.Div(bOvO + rvbGvG);
else
    oil = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 );
    divOil = s.Div(bOvO);
end

% Conservation of mass for gas. Again, we have two cases depending on
% whether the model allows us to dissolve the gas phase into the oil phase.
if model.disgas
    % The gas transported in the oil phase.
    rsbOvO = s.faceUpstr(upco, rs).*bOvO;
    
    gas = (s.pv/dt).*( pvMult.* (bG.* sG  + rs.* bO.* sO) - ...
                       pvMult0.*(bG0.*sG0 + rs0.*bO0.*sO0 ));
    divGas = s.Div(bGvG + rsbOvO);
else
    gas = (s.pv/dt).*( pvMult.*bG.*sG - pvMult0.*bG0.*sG0 );
    divGas = s.Div(bGvG);
end

% Put the set of equations into cell arrays along with their names/types.
if model.water
    eqs      = {water, oil, gas};
    divTerms = {divWater, divOil, divGas};
    names    = {'water', 'oil', 'gas'};
    types    = {'cell', 'cell', 'cell'};
    rho      = {rhoW, rhoO, rhoG};
    mob      = {mobW, mobO, mobG};
    pressures = {pW, p, pG};
else
    eqs      = {oil, gas};
    divTerms = {divOil, divGas};
    names    = {'oil', 'gas'};
    types    = {'cell', 'cell'};
    rho      = {rhoO, rhoG};
    mob      = {mobO, mobG};
    pressures = {p, pG};
end
% Add in any fluxes / source terms prescribed as boundary conditions.
dissolved = model.getDissolutionMatrix(rs, rv);

[eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                 pressures, sat, mob, rho, ...
                                                 dissolved, {}, ...
                                                 drivingForces);

% Add in and setup well equations
[eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, types, wellSol0, wellSol, wellVars, wellMap, p, mob, rho, dissolved, {}, dt, opt);

% Finally, adding divergence terms to equations
for i = 1:numel(divTerms)
    eqs{i} = eqs{i} + divTerms{i};
end

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end

