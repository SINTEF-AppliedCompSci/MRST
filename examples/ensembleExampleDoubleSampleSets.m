%% Ensemble of 3D reservoirs
% Whereas this also serves as a minimal example, it demonstrates how to use
% different sample classes and combine them so that an ensemble has
% uncertain parameters related to different properties in the model
% (uncertain rock and uncertain well properties).
% To show this we create a small ensemble of a 3D reservoir modelling two
% phase flow (oil and water) between two injectors and two producers.
mrstModule add ad-core ad-blackoil mrst-gui ad-props ...
    example-suite incomp ensemble 

mrstVerbose off

%% Choose a base problem and ensemble size
% The base problem contains all properties that are common throughout the
% ensemble, and here we have implemented our own super simple example (see 
% example_template.m for the MRSTExample template).

baseProblemName = 'ensemble_base_problem_3d_reservoir';
baseProblemOptions = {};


% Change these flags to investigate the baseExample
simulateExample = false;
plotSimulation = false;
rerunBaseProblemFromScratch = false;

baseExample = MRSTExample(baseProblemName, baseProblemOptions{:});

if simulateExample
    problem = baseExample.getPackedSimulationProblem();
    if rerunBaseProblemFromScratch
        clearPackedSimulatorOutput(problem);
    end
    simulatePackedProblem(problem);

    [wellSols, states, reports] = getPackedSimulatorOutput(problem);
    if plotSimulation
        baseExample.plot(states);
    end
end

%% an ensemble of stochastic rock realizations 
% RockSamples is a superclass of BaseSamples, and thereby implements all
% the sample functionality required to map new rock data to the base
% problem. Here, we precompute 20 rock realizations and use these as
% stochastic samples.

ensembleSize = 20;

rockData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    rockData{i}.poro = gaussianField(baseExample.model.G.cartDims, [0.2 0.4]); 
    rockData{i}.perm = rockData{i}.poro.^3.*(1e-5)^2./(0.81*72*(1-rockData{i}.poro).^2);
end

rockSamples = RockSamples('data', rockData);

%% Select quantity of interest class
% We here demonstrates how the WellQoI class can be used to store
% production information for two fields for both producers.

qoi = WellQoI(...
    'wellNames', {'P1', 'P2'}, ...
    'fldname', {'qOs', 'qWs'}, ...
    'cumulative', false);


%% Create ensemble
% Combining the baseExample, the rock samples and the well QoI.
% We will run the ensemble in parallel across 8 workers, and we choose to
% delete any existing simulation results and run the ensemble from scratch

rockEnsemble = MRSTEnsemble(baseExample, rockSamples, qoi, ...
    'simulationStrategy', 'background', ...
    'maxWorkers', 8, ...
    'reset', true, ...
    'verbose', true);

%% Run ensemble
rockEnsemble.simulateEnsembleMembers('plotProgress', true);

%% Plot results
% Set optional input argument `subplots` to `true` to plot each QoI field
% for all wells in the same figure on top of each other (`subplotDir =
% 'vertical'`) or side-by-side ((`subplotDir =  'horizontal'`)
close all, rockEnsemble.plotQoI('subplots', true, 'subplotDir', 'vertical');

%% Create another ensemble using stochastic well indices.
% The class WellSamples are used in the same way as RockSamples, and are
% also a superclass of BaseSamples. Instead of rock properties, however,
% its data property now holds well properties. In this example, we have
% stochastic well production indices (WI).

wellSampleData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    wellSampleData{i}.WI = rand(1,4)*1e-11;
end

wellSamples = WellSamples('data', wellSampleData);

%% Define new ensemble
wellEnsemble = MRSTEnsemble(baseExample, wellSamples, qoi, ...
    'simulationStrategy', 'background', ...
    'maxWorkers', 8, ...
    'reset', true, ...
    'verbose', true);

%% Simulate and plot
wellEnsemble.simulateEnsembleMembers('plotProgress', true);

%% Plot result
close all, wellEnsemble.plotQoI('subplots', true, 'subplotDir', 'vertical');

%% Create another ensemble that combines both sampling strategies
% Finally, we want to run an ensemble that has stochastic rock properties
% AND stochastic well properties. To achieve this, we have implemented the
% class WellRockSamples that is a superclass of both WellSamples and
% RockSamples. 
% The data property in this class expects to have a cell array of structs
% with fields rock and well, both these being structs with rock properties
% and well properties, respectively.
comboData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    comboData{i}.rock = rockData{i};
    comboData{i}.well = wellSampleData{i};
end
comboSamples = WellRockSamples('data', comboData);

%% Create another ensemble that combines both sampling strategies
% Finally, we want to run an ensemble that has stochastic rock properties
% AND stochastic well properties. This is conveniently implemented in the
% `CompositeSamples` class, where we can combine any number of parent
% samples. The optional input argument `tensorProduct = true` uses all
% tensor products of all parent samples. In effect, this means that we use a
% different seed to get each parent sample. If this is false, we use the
% same seed to get all parent samples.
compSamples = CompositeSamples({rockSamples, wellSamples}, 'tensorProduct', false);

%% Define new ensemble
compEnsemble = MRSTEnsemble(baseExample, compSamples, qoi, ...
    'simulationStrategy', 'background', ...
    'maxWorkers', 8, ...
    'reset', true, ...
    'verbose', true);

%% Simulate and plot
compEnsemble.simulateEnsembleMembers('plotProgress', true);

%% Plot results
color = lines(2); color = color(2,:);
close all, compEnsemble.plotQoI('subplots'  , true      , ...
                                'subplotDir', 'vertical', ...
                                'color'     , color     );

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.
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
