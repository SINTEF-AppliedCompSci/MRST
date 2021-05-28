%% 2D Three-Phase Surfactant-Polymer Injection Case
% Test of the generic and legacy solvers for EOR
mrstModule add ad-core ad-blackoil ad-eor ad-props ...
               deckformat mrst-gui

%% Set up model and initial conditions
% The data required for the example
% The following are all the files needed for this tutorial
% Two files are the data for the simulation of surfactant polymer flooding.
fn = fullfile(mrstPath('ad-eor'), 'examples', 'surfactant-polymer', 'SURFACTANTPOLYMER2D.DATA');
% Construct physical model, initial state and dynamic well controls.
[state0, model, schedule, nls] = ...
   initEclipseProblemAD(fn, 'timestepstrategy', 'none');
arg = {model.G, model.rock, model.fluid, ...
      'disgas', model.disgas, 'vapoil', model.vapoil,...
      'inputdata', model.inputdata};
%%  Run with three-phase surfactant polymer model
lmodelsp = ThreePhaseSurfactantPolymerModel(arg{:});
scheduleSP = schedule;
[wellSolsSP, statesSP, reportsSP] = ...
    simulateScheduleAD(state0, lmodelsp, scheduleSP, ...
                    'NonLinearSolver', nls);

%% Run with generic surfactant-polymer model
% Model can use any combination of phases and components.
gmodelsp = GenericSurfactantPolymerModel(arg{:});
scheduleGSP = schedule;

[wellSolsGSP, statesGSP, reportsGSP] = ...
    simulateScheduleAD(state0, gmodelsp, scheduleGSP, ...
                    'NonLinearSolver', nls);

%% Run just the black-oil case with no EOR effects considered
bomodel = GenericBlackOilModel(arg{:});
[wellSolsBO, statesBO, reportsBO] = ...
    simulateScheduleAD(state0, bomodel, schedule, ...
                    'NonLinearSolver', nls);

%% Compare results
plotWellSols({wellSolsGSP, wellSolsSP, wellSolsBO}, 'datasetnames', {'GenericSP', 'OldModel', 'Blackoil'});
%% Plot state function diagrams
figure;
plotStateFunctionGroupings(bomodel)
title('Black-oil')
figure;
plotStateFunctionGroupings(gmodelsp)
title('EOR')

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
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
