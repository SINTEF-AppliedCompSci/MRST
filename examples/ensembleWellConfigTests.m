%% Ensemble of 3D reservoirs
% This is a test to see that different types of well configurations work
% with ensembles.
% The example uses a small 3D reservoir with 4 wells. The wells can
% perforate either a single cell or several cells, and the number of
% perforated cells can be different between each well). We also ensure that
% we can handle a schedule that varies with time.
mrstModule add ad-core ad-blackoil mrst-gui ad-props ...
    test-suite incomp ensemble 

mrstVerbose off

%% Choose a base problem and ensemble size
% The base problem contains all properties that are common throughout the
% ensemble, and here we have implemented our own super simple example (see 
% example_template.m for the MRSTExample template).

baseProblemName = 'ensemble_base_problem_3d_reservoir';
standardProblemOptions = {};
longWellsProblemOptions = {'longWells', true, ...
                           'randomSchedule', false};
randomScheduleProblemOptions = {'longWells', false, ...
                                'randomSchedule', true};

% Change these flags to investigate the baseExample
simulateExample = false;
plotSimulation = true;
rerunBaseProblemFromScratch = true;

baseExample = TestCase(baseProblemName, randomScheduleProblemOptions{:});

if simulateExample
    problem = baseExample.getPackedSimulationProblem();
    if rerunBaseProblemFromScratch
        clearPackedSimulatorOutput(problem);
    end
    simulatePackedProblem(problem);

    [wellSols, states, reports] = getPackedSimulatorOutput(problem);
    if plotSimulation
        baseExample.plot(states);
        plotWellSols(wellSols);
    end
end

%% Select quantity of interest class
% We here demonstrates how the WellQoI class can be used to store
% production information for two fields for both producers.

qoi = WellQoI(...
    'wellNames', {'P1', 'P2'}, ...
    'names', {'qOs', 'qWs'}, ...
    'cumulative', false);


%% Create ensemble using stochastic well indices.
% The class WellSamples are used in the same way as RockSamples, and are
% also a superclass of BaseSamples. Instead of rock properties, however,
% its data property now holds well properties. In this example, we have
% stochastic well production indices (WI).

ensembleSize = 20;
wellSampleData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    wellSampleData{i}.WI = rand(1,4)*1e-11;
end

wellSamples = WellSamples('data', wellSampleData);

%% Define new ensemble
wellEnsemble = MRSTEnsemble(baseExample, wellSamples, qoi, ...
    'simulationStrategy', 'serial', ...
    'maxWorkers', 8, ...
    'reset', true, ...
    'verbose', true);

%% Simulate and plot
close all;
wellEnsemble.simulateEnsembleMembers('subplots', true,....
                                     'subplotDir', 'vertical',...
                                     'clearFigure', true);

%% Plot result
close all, wellEnsemble.plotQoI('subplots', true, 'subplotDir', 'vertical');


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
