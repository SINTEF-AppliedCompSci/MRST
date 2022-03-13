%% Norne: test calibrated GPSNet model versus different schedules
% The purpose of this script is to assess to what extent a GPSNet trained
% as described in the script norneGPSNetAdjoint generalizes to simulations
% run with well controls that are different from the ones used to train the
% network model. The parameter caseNo determines how the well controls are
% setup:
%       0 - no modification of the original prediction case
%       1 - 25% random perturbation of well rates and 10% perturbation of
%           bhp controls
%       2 - shut in the dominant producer (P1) after 2/3 of the simulation
%           time (8 years)
%       3 - convert injectors I1/I2 to producers and shut in P1
%       4 - shut in P1 during the middle third of the simulation horizon
%           (i.e., in the period from year 4 to year 8)
%
% To save you from having to train the GPSNet over many iterations (which
% can take considerable time), we have stored parameters after 50, 100,
% 200, 300, 400, 500, and 750 iterations for GPSNet with all-to-all
% connections and injector-to-producer settings. Notice, however, that
% these have been trained for a particular petrophysical realization of
% Norne may not work well on different realizations. 
mrstModule add ad-core ad-blackoil deckformat diagnostics...
               mrst-gui ad-props incomp optimization...
               network-models test-suite linearsolvers 

if ~exist('caseNo','var'), caseNo = 1; end

%% Setup 3D reference model
predCase  = TestCase('norne_simple_wo');
predCase  = modifyNorneTest(predCase, caseNo);
predProbl = predCase.getPackedSimulationProblem();
%clearPackedSimulatorOutput(trueProb)
simulatePackedProblem(predProbl);

[wellSolPred, statesPred] = getPackedSimulatorOutput(predProbl);
predModel     = predCase.model;
predSchedule  = predProbl.SimulatorSetup.schedule;
Wpred         = predSchedule.control.W;

%% Create the GPSNet model
% First, we create the network
Wnw = Wpred;
for i = 1:numel(Wnw)
    num_cells = numel(Wnw(i).cells);
    Wnw(i).cells = Wnw(i).cells(round(num_cells/2));
end
%ntwrk =  Network(Wnw, predModel.G, 'type', 'injectors_to_producers', ...
%                 'injectors', 1:6, 'producers', 7:11);
ntwrk =  Network(Wnw, predModel.G, 'type', 'all_to_all');

gravity off
gpsNet    = GPSNet(predModel, ntwrk, Wpred, 10);
predSetup = gpsNetSimulationSetup(gpsNet, predSchedule);

%% Specify parameters and define mismatch
[cellEdgeNo, cellIx] = gpsNet.getMapping('cells');
[faceEdgeNo, faceIx] = gpsNet.getMapping('faces');
config = {
    ...%name      include     scaling    boxlims     lumping     subset   relativeLimits
    'porevolume',       1,   'linear',       [],  cellEdgeNo,   cellIx,   [.001 4]
    'conntrans',        1,   'log',          [],          [],       [],   [.001 100]
    'transmissibility', 1,   'log'   ,       [],  faceEdgeNo,   faceIx,   [.001 100]};
params = [];
for k = 1:size(config,1)
    if config{k, 2} == 0, continue, end
    params = addParameter(params, predSetup, ...
        'name',    config{k,1}, 'scaling', config{k,3}, ...
        'boxLims', config{k,4}, 'lumping', config{k,5}, ...
        'subset',  config{k,6}, 'relativeLimits',config{k,7});
end

weighting  = {'WaterRateWeight',  day/10000, ...
              'OilRateWeight',    day/20000, ...
              'BHPWeight',        1/(500*barsa)};
mismatchFn = @(model, states, schedule, states_ref, tt, tstep, state) ...
    matchObservedOW(model, states, schedule, states_ref,...
                   'computePartials', tt, 'tstep', tstep, weighting{:},...
                   'state',state,'from_states',false);
               
%% Load parameters and evaluate mismatch
% We run the GPSNet with parameters obtained after an increasing number of
% iterations. These are stored in a file, which contains two different
% parameter sets: 
%  poptAll     - GPSNet set up with all-to-all connections
%  poptInjProd - GPSNet set up with injector-to-producer connections
its = [50 100:100:500 750]; misfit = 0*its;
data = load('gpsnet-norne.mat');
for n=1:numel(its) 
    [misfit(n),~,wellSol] = evaluateMatch(data.poptAll(:,n), ...
        mismatchFn, predSetup, params, statesPred, 'Gradient', 'none');
end
for n=1:numel(its), fprintf('%d\t%.2e\n', its(n), -misfit(n)); end

%% Plot well curves
plotWellSols({wellSolPred,wellSol}, ...
    {predSchedule.step.val,predSetup.schedule.step.val}, ...
    'datasetnames', {'reference','GPSNet'}, 'zoom', true, ...
    'field', 'qOs', 'SelectedWells', 7);
set(gcf,'Name',predCase.name);

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
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
