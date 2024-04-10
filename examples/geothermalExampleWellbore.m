%% Geothermal wellbore model example
% Example demonstrating how to use WellboreModel and discrete fracture
% modelling (DFM) in geothermal simulations

%% Add modules
mrstModule add test-suite
mrstModule add ad-core ad-props
mrstModule add geothermal compositional dfm
mrstModule add wellpaths upr
mrstModule add mrst-gui
mrstModule add linearsolvers
mrstVerbose on

%% Set up test case
% We use a simple fivespot pattern with horizontal fractures
cartDims = [11,11,24];
test     = TestCase('fivespot_geothermal', 'cartDims', cartDims, 'dfm', true);

% Remove facility model and set up groups
test.model.FacilityModel = [];
groups = [];
hot  = {'Hot'};
cold = {'Cold-1-1', 'Cold-2-1', 'Cold-2-2', 'Cold-1-2'};
groups = addFacilityGroup(groups, hot , 'name', 'HotGroup');
groups = addFacilityGroup(groups, cold, 'name', 'ColdGroup');
test = convertToReservoirWellboreModel(test, 'groups', groups);

% Define controls on group level
% Charging
ctrl = test.schedule.control(1).Wellbore;
[ctrl.W.type] = deal('group');
ctrl.groups = groups;
ctrl.groups(1).type = 'rate';
ctrl.groups(2).type = 'bhp';
val = sum(vertcat(ctrl.W(5).val));
ctrl.groups(1).val = val;
ctrl.groups(2).val = ctrl.W(1).val;
ctrl.groups(1).T = ctrl.W(5).T;
ctrl.groups(2).T = ctrl.W(1).T;
test.schedule.control(1).Wellbore = ctrl;
% Discharging
ctrl = test.schedule.control(2).Wellbore;
[ctrl.W.type] = deal('group');
ctrl.groups = groups;
ctrl.groups(1).type = 'bhp';
ctrl.groups(2).type = 'rate';
val = vertcat(ctrl.W(1:4).val);
ctrl.groups(1).val = ctrl.W(5).val;
ctrl.groups(2).val = sum(val);
ctrl.groups(1).T = ctrl.W(5).T;
ctrl.groups(2).T = ctrl.W(1).T;
test.schedule.control(2).Wellbore = ctrl;

%% Set up linear solver
lsol = setUpGeothermalWellboreLinearSolver(test.model);
test.model.G.cells.num = test.model.submodels.Reservoir.G.cells.num + ...
nnz(test.model.submodels.Wellbore.G.cells.type == 0 & ...
    test.model.submodels.Wellbore.G.cells.hybrid == 0);

%% Simulate
problem = test.getPackedSimulationProblem('LinearSolver', lsol);
simulatePackedProblem(problem, 'restartStep', 1);

%% Simulate reference case without wellbore model
testRef = TestCase('fivespot_geothermal', 'cartDims', cartDims, 'dfm', true);
problemRef = testRef.getPackedSimulationProblem();
simulatePackedProblem(problemRef, 'restartStep', 1);

%% Compare
% Get simulation output
[~, states, reports] = getPackedSimulatorOutput(problem);
% Extract wellSols from state
wellSols = getWellSolsFromWellboreState(test.model, states);
[wellSolsRef, statesRef, reportsRef] = getPackedSimulatorOutput(problemRef);

% Plot simulations with wellbore and reference
test.plot(states);
testRef.plot(statesRef);
% Plot difference
dstates = cellfun(@(st1, st2) ...
    compareStructs(st1.Reservoir, st2, ...
    'includeStructs', false, 'relative', true, 'verbose', false), ...
    states, statesRef, 'UniformOutput', true ...
);
test.plot(dstates);

% Plot well solutions
plotWellSols({wellSols, wellSolsRef}, test.schedule.step.val);

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>