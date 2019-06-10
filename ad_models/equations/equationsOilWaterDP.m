function [problem, state] = equationsOilWaterDP(state0, state, model, dt, drivingForces, varargin)
% Generate linearized problem for the two-phase oil-water model
%
% SYNOPSIS:
%   [problem, state] = equationsOilWater(state0, state, model, dt, drivingForces)
%
% DESCRIPTION:
%   This is the core function of the two-phase oil-water solver. This
%   function assembles the residual equations for the conservation of water
%   and oil as well as required well equations. By default, Jacobians are
%   also provided by the use of automatic differentiation.
%
% REQUIRED PARAMETERS:
%   state0    - Reservoir state at the previous timestep. Assumed to have
%               physically reasonable values.
%
%   state     - State at the current nonlinear iteration. The values do not
%               need to be physically reasonable.
%
%   model     - TwoPhaseOilWaterModel-derived class. Typically,
%               equationsOilWater will be called from the class
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
%   problem - LinearizedProblemAD class instance, containing the water
%             and oil conservation equations, as well as well equations
%             specified by the WellModel class.
%
%   state   - Updated state. Primarily returned to handle changing well
%             controls from the well model.
%
% SEE ALSO:
%   equationsBlackOil, TwoPhaseOilWaterModel

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
             'reverseMode', false,...
             'resOnly', false,...
             'iteration', -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.W;
s = model.operators;

%% Dual porosity: adding matrix
% Properties at current timestep
[p, sW, wellSol, pm, swm] = model.getProps(state, 'pressure', 'water', 'wellsol',...
                                         'pressure_matrix', 'water_matrix');
% Properties at previous timestep
[p0, sW0, wellSol0, pm0, swm0] = model.getProps(state0, 'pressure', 'water', 'wellSol',...
                                             'pressure_matrix', 'water_matrix');

%%
[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

% Initialize independent variables.
if ~opt.resOnly
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode
        [p, sW, pm, swm, wellVars{:}] = model.AutoDiffBackend.initVariablesAD(p, sW, pm, swm, wellVars{:});
    else
        wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
        [p0, sW0, pm0, swm0, wellVars0{:}] = model.AutoDiffBackend.initVariablesAD(p0, sW0, pm0, swm0, wellVars0{:}); %#ok
    end
end
primaryVars = {'pressure', 'sW', 'pressure_matrix', 'swm', wellVarNames{:}};
% We will solve for pressure, water saturation (oil saturation follows via
% the definition of saturations) and well rates + bhp.

%% Fracture
sO  = 1 - sW;
sO0 = 1 - sW0;

sat = {sW, sO};
sat0 = {sW0, sO0};

% Update state with AD-variables
state = model.setProps(state, {'s', 'pressure'}, {sat, p});
state0 = model.setProps(state0, {'s', 'pressure'}, {sat0, p0});
% Set up properties
state = model.initStateFunctionContainers(state);

[b, pv] = model.getProps(state, 'ShrinkageFactors', 'PoreVolume');
[b0, pv0] = model.getProps(state0, 'ShrinkageFactors', 'PoreVolume');
[phaseFlux, flags] = model.getProps(state, 'PhaseFlux',  'PhaseUpwindFlag');

[bW, bO] = deal(b{:});
[bW0, bO0] = deal(b0{:});
[vW, vO] = deal(phaseFlux{:});
[upcw, upco] = deal(flags{:});

[pressures, mob, rho] = model.getProps(state, 'PhasePressures', 'Mobility', 'Density');

%% Matrix
fm = model.fluid_matrix;
%check for p-dependent porv mult:
pvMultm = 1; pvMultm0 = 1;
if isfield(fm, 'pvMultR')
    pvMultm =  fm.pvMultR(pm);
    pvMultm0 = fm.pvMultR(pm0);
end
pvm = pvMultm.*model.operators.pv_matrix;
pvm0 = pvMultm0.*model.operators.pv_matrix;

pcOWm = 0;
pcOWm0 = 0;
if isfield(model.fluid_matrix, 'pcOW') && ~isempty(swm)
    pcOWm  = model.fluid_matrix.pcOW(swm);
    pcOWm0  = model.fluid_matrix.pcOW(swm0);
end
pwm = pm - pcOWm;
pwm0 = pm0 - pcOWm0;
bWm = fm.bW(pwm);
bWm0 = fm.bW(pwm0);
bOm = fm.bO(pm);
bOm0 = fm.bO(pm0);

%% Transfer
vb = model.G.cells.volumes;

matrix_fields.pom = pm;
matrix_fields.swm = swm;
fracture_fields.pof = p;
fracture_fields.swf = sW;

transfer_model = model.transfer_model_object;

[Talpha] = transfer_model.calculate_transfer(model,fracture_fields,matrix_fields);

Twm = vb.*Talpha{1};
Tom = vb.*Talpha{2};
%% Output additional information
if model.outputFluxes
    state = model.storeFluxes(state, vW, vO, []);
end
if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, []);
    state = model.storeMobilities(state, mob{:}, []);
    state = model.storeUpstreamIndices(state, upcw, upco, []);
end

% EQUATIONS ---------------------------------------------------------------
% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.
%% Fracture
bOvO = s.faceUpstr(upco, bO).*vO;
bWvW = s.faceUpstr(upcw, bW).*vW;

% Conservation of mass for water
water_fracture = (1/dt).*( pv.*bW.*sW - pv0.*bW0.*sW0 )  + Twm;
% Conservation of mass for oil
oil_fracture = (1/dt).*( pv.*bO.*sO - pv0.*bO0.*sO0 )  + Tom;

%% Matrix
% Conservation of mass for water
water_matrix = (1/dt).*( pvm.*bWm.*swm - pvm0.*bWm0.*swm0 ) - Twm;
% Conservation of mass for oil
oil_matrix = (1/dt).*( pvm.*bOm.*(1-swm) - pvm0.*bOm0.*(1-swm0) ) - Tom;


%%
eqs = {water_fracture, oil_fracture, water_matrix, oil_matrix};
names = {'water', 'oil', 'water_matrix', 'oil_matrix'};
types = {'cell', 'cell', 'cell', 'cell'};

% Add in any fluxes / source terms prescribed as boundary conditions.
[eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                                 pressures, sat, mob, rho, ...
                                                                 {}, {}, ...
                                                                 drivingForces);
% Finally, add in and setup well equations
[eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, types, wellSol0, wellSol, wellVars, wellMap, p, mob, rho, {}, {}, dt, opt);
% Add in fluxes
eqs{1} = eqs{1} + s.Div(bWvW);
eqs{2} = eqs{2} + s.Div(bOvO);


problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end

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
