function [problem, state] = equationsCO2(state0, state, model, dt,...
                                                   drivingForces, varargin)
% Function to assemble the linearized equations for the co2-water system.
%
% This function is modified from a file in The MATLAB Reservoir Simulation
% Toolbox (MRST), see
%   mrst/modules/ad-eor/utils/equationsOilWaterPolymer.m
%
% We refer to that function for a complete commented version of the file.
% In this file we comment on some of the lines. 

%{
Partial copyright 2009-2020, SINTEF Digital, Mathematics & Cybernetics.
Partial copyright 2020, NORCE Norwegian Research Centre AS, Computational
Geosciences and Modeling.

This file is part of the ad-micp module.

ad-micp is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ad-micp is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%}
opt = struct('Verbose', mrstVerbose, 'reverseMode', false, 'resOnly',...
                                                   false, 'iteration', -1);
opt = merge_options(opt, varargin{:});
op = model.operators;
rhoW = model.fluid.rhoWS*model.fluid.cells;
rhoO = model.fluid.rhoOS*model.fluid.cells;

% Properties at current timestep
[p, sW, wellSol] = model.getProps(state, 'pressure', 'water', 'wellSol');
% Properties at previous timestep
[sW0, wellSol0] = model.getProps(state0, 'water', 'wellSol');

[wellVars, wellVarNames, wellMap] = ...
                       model.FacilityModel.getAllPrimaryVariables(wellSol);
% Initialize independent variables.
if ~opt.resOnly
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode
        [p, sW, wellVars{:}] = initVariablesADI(p, sW, wellVars{:});
        primaryVars = [{'pressure'}, {'sW'}, wellVarNames(:)'];
    else
        [sW0] = initVariablesADI(sW0);
        primaryVars = {'pressure', 'sW'};
    end
else
    primaryVars = {'pressure', 'sW'};
end

%sO -> co2 saturation
sO  = 1 - sW;
sO0 = 1 - sW0;

% Evaluate relative permeabilities
krW = sW;
krO = sO;

T = op.T;
gdz = model.getGravityGradient();

[vW, mobW] = getFluxAndPropsWater(model, p, krW, T, gdz);
[vO, mobO] = getFluxAndPropsCO2(model, p, krO, T, gdz);

if model.outputFluxes
    state = model.storeFluxes(state, vW, vO, 0*vW);
end

% Conservation of mass for water
water = (op.pv/dt).*(sW - sW0) + op.Div(vW);
% Conservation of mass for co2
oil = (op.pv/dt).*(sO - sO0) + op.Div(vO);

eqs   = {water, oil};
names = {'water', 'oil'};
types = {'cell', 'cell'};

rho = {rhoW, rhoO};
mob = {mobW, mobO};
sat = {sW, sO};

[eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types,...
                     state, {p, p}, sat, mob, rho, {}, {}, drivingForces);

[eqs, names, types, state.wellSol] = model.insertWellEquations(eqs,...
                              names, types, wellSol0, wellSol, wellVars,...
                                    wellMap, p, mob, rho, {}, {}, dt, opt);

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end
