function [problem, state] = equationsOilWater(state0, state, model, dt, drivingForces, varargin)
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
opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'iteration', -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.W;
s = model.operators;

% Properties at current timestep
[p, sW, wellSol] = model.getProps(state, 'pressure', 'water', 'wellsol');
% Properties at previous timestep
[p0, sW0, wellSol0] = model.getProps(state0, 'pressure', 'water', 'wellSol');

[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

% Initialize independent variables.
if ~opt.resOnly
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode
        [p, sW, wellVars{:}] = model.AutoDiffBackend.initVariablesAD(p, sW, wellVars{:});
    else
        wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
        [p0, sW0, wellVars0{:}] = model.AutoDiffBackend.initVariablesAD(p0, sW0, wellVars0{:}); %#ok
    end
end
primaryVars = {'pressure', 'sW', wellVarNames{:}};
% We will solve for pressure, water saturation (oil saturation follows via
% the definition of saturations) and well rates + bhp.
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
bOvO = s.faceUpstr(upco, bO).*vO;
bWvW = s.faceUpstr(upcw, bW).*vW;

% Conservation of mass for water
water = (1/dt).*( pv.*bW.*sW - pv0.*bW0.*sW0 );

% Conservation of mass for oil
oil = (1/dt).*( pv.*bO.*sO - pv0.*bO0.*sO0 );

eqs = {water, oil};
names = {'water', 'oil'};
types = {'cell', 'cell'};

% Add in any fluxes / source terms prescribed as boundary conditions.
[eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                                 pressures, sat, mob, rho, ...
                                                                 {}, {}, ...
                                                                 drivingForces);

% Add aquifer contributions if any.
if ~isempty(model.AquiferModel)
    eqs = addAquifersContribution(model.AquiferModel, eqs, names, state, dt);
end

% Finally, add in and setup well equations
[eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, types, wellSol0, wellSol, wellVars, wellMap, p, mob, rho, {}, {}, dt, opt);
% Add in fluxes
eqs{1} = eqs{1} + s.Div(bWvW);
eqs{2} = eqs{2} + s.Div(bOvO);


problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end

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
