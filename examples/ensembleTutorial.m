%% Ensemble tutorial
% This tutorial goes thorugh how to set up and run an ensemble simulation
% in MRST using the ensemble module

%% Add necessary modules
mrstModule add ad-core ad-props ad-blackoil 
mrstModule add ensemble
mrstModule add test-suite
mrstModule add mrst-gui
mrstVerbose on

%% Set up base case
baseCase = TestCase('qfs_wo');
problem = baseCase.getPackedSimulationProblem();

%% Generate an ensemble
% The stochastic component (or uncertain parameter) of any ensemble
% simulation is implemented in a sample class (a superclass of 
% `BaseSample`), and a specific realization of such a stochastic component
% is referred to as a sample. We can set up a sample class in three
% different ways, by providing
% 1) a function that generate stochastic samples;
% 2) a cell array of precomputed samples;
% 3) an instance of ResultHandler that points to a location where
%    precomputed samples are stored.
% We will show you all three possibilities in this example.

%% Using a generator function
% We will use the function 'generateRockSample', which generates rock
% samples based on a stationary Gaussian process on the [0,1]^d box [1].
% Inputs to the generator function in the sample class should be the
% problem setup and seed for controlling the random number generator, so we
% make a function handle taking in these arguments.
generatorFn = @(problem, seed) ...
    generateRockSample(problem.SimulatorSetup.model.G.cartDims, ...
                                           'seed', seed, 'toVector', true);

% The class RockSamples implements routines for getting stochastic rock
% realizations, and setting such a sample to the model. The latter also 
% includes updating all model operators depending on the rock.
samplesFn = RockSamples('generatorFn', generatorFn);

% Notice that the data property of samplesFn is empty, the
% generatorFn property in non-empy, whereas the num propery says inf. The
% latter refers to the number of samples, and is inf since we in principle
% can generate as many samples as we want.
disp(samplesFn);

%% Using a cell array of precomputed samples
% To illustrate how we can use precomputed rock samples, we use the same
% function to generate an ensemble of 50 realizations.
ensembleSize = 50;
data         = cell(ensembleSize, 1);
for i = 1:ensembleSize
    data{i} = generatorFn(problem, i);
end
samplesCell = RockSamples('data', data);

% This time, the 'generatorFn' property is empty, whereas the 'data'
% property is a cell array of the 50 realizations we just made. Note also
% that the 'num' property says 50.
disp(samplesCell);

%% Using a ResultHandler
% The problem size is very often so large that holding all ensemble members
% in memory is intractable. In cases when we have a set of samples that are
% stored to file, we can provide a ResultHandler to the RockSample class
% that facilitates loading the samples from disk.
dataDir = fullfile(mrstOutputDirectory(), 'ensemble', 'tutorial');
if ~exist(dataDir, 'dir'), mkdir(dataDir); end
rh = ResultHandler('dataDirectory', dataDir  , ... % Example root folder
                   'dataFolder'   , 'samples', ... % Sample folder
                   'dataPrefix'   , 'sample_');    % samples_<seed>.mat
rh(1:ensembleSize) = data;           % Store samples to file
samplesRH = RockSamples('data', rh); % Set up rock samples using the result handler

% As in the previous case, the 'generatorFn' property is empty and the
% 'num' property says 50, but the 'data' property is now the ResultHandler
% we just created.
disp(samplesRH);

%% Inspect ensemble subset
% Before creating the complete ensemble, we investigate a subset of the
% samples define above, and illustrate how we obtain a single realization
% of the samples in terms of a prolem.
data    = samplesCell.getSample(13, problem);   % Get sample numer 13
problem = samplesCell.setSample(data, problem);  % Set sample to problem
% Like the rock structure of a model, the sample has a perm and poro field
disp(data);

% Inspect rock samples
baseCase.plot(problem.SimulatorSetup.model.rock, 'log10', true); colormap(pink);
for i = 1:5
    data    = samplesCell.getSample(i, problem);   % Get sample number i
    problem = samplesCell.setSample(data, problem); % Set sample to problem
    % Inspect rock sample
    baseCase.plot(problem.SimulatorSetup.model.rock, 'log10', true); colormap(pink);
end

%% Quantity of interest
% The quantity that we are interested in from an ensemble simulation is
% called quantity of interest (QoI). In this example, we will look at two
% QoIs: the water saturation after 20, 300, and 700 days of injection, and
% the oil production rate over the entire simulation.

%% Reservoir state QoIs
% All QoIs derived from the reservoir state (except wells) are implemented
% in ReservoirStateQoI. The 'name' input refers to the name of a
% property that can be fetched or computed from the state using
% model.getProp, whereas the 'time' input is the time at which we want the
% quantity. This may be a vector, and for each element, the class will pick
% the timestep correponding to the closest timestamp.
qoiState = ReservoirStateQoI('names', 'sW', 'time', [20, 300, 700]*day);
disp(qoiState);

%% Well output QoIs
% All QoIs related to well output are implemented in WellQoI. This class
% uses the function getWellOutput to get the results from the well
% solutions, and therefore has properties inputs of this function:
% 'fldnames' is field name of interest stored on wellSol, whereas 'wellIndices'
% is a vector with the well number(s). An alternative to 'wellIndices' is to 
% provide 'wellNames' instead as a cell array of well names.
qoiWell = WellQoI('names', 'qOs', 'wellIndices', 2);
disp(qoiWell);

%% Running a single sample
% We have now defined all ingredients necessary to set up an instance of
% the ensemble class. First, however, we set up and simulate a single
% sample to illustrate the steps required to obtain a single problem
% realization, running the problem, and obtaining the quantity of interest.

problem = baseCase.getPackedSimulationProblem(); % Get packed problem
data    = samplesCell.getSample(13, problem);   % Get sample numer 13
% Like the rock structure of a model, the sample has a perm and poro field
disp(data);
problem = samplesCell.setSample(data, problem);  % Set sample to problem
% Inspect rock sample
baseCase.plot(problem.SimulatorSetup.model.rock, 'log10', true); colormap(pink);
simulatePackedProblem(problem); % Simulate

% To inspect the quantities of interest directly, we first need to
% match the configurations of the QoI objects with the configurations of
% the problem at hand. Note that this step is done automatically within the 
% MRSTEnsemble class constructor.
qoiStateValidated = qoiState.validateQoI(problem);
qoiWellValidated  = qoiWell.validateQoI(problem);
disp(qoiStateValidated);
disp(qoiWellValidated);

% Quantities of interest are computed from the problem
stateData = qoiStateValidated.computeQoI(problem); % Water saturation after 20, 300, and 700 days
wellData  = qoiWellValidated.computeQoI(problem);  % Producer oil production rate

%% Plot the results
% The results can be plotted using the plotting functionality found in the
% example, or simply by creating a plot from scratch.
close all
% Water saturation
baseCase.plot(stateData); colormap(bone);

% Oil production rate
time = cumsum(baseCase.schedule.step.val)/day;
figure(), plot(time, wellData.qOs*day, 'LineWidth', 1.5); % Convert to days
xlim([0, time(end)]), box on, grid on, xlabel('Time (days)');

%% Plot the results using QoI functionality
% ... or, we can plot the results using the plotting functions in the QoI
% classes.
close all
qoiWellValidated.plotQoI(baseCase, wellData);
baseCase.figure();
qoiStateValidated.plotQoI(baseCase, stateData);

%% Set up the ensemble
% The MRSTEnsemble class conveniently combines a problem setup, a
% samples object and a QoI, and implements the functionality for everything
% related to setting up and simulating the ensemble members and computing 
% their corresponding QoI. The first input parameter can either be an
% instance of an TestCase, or the name of an TestCase function. In 
% the latter case, MRSTEnsemble will first set up the example, and any 
% input arguments not used by the MRSTEnsemble class itself will be passed
% on to the MRSTExample class when creating the example instance.
%
% The ensemble can run simulations in three different ways, depending on
% the optional input parameter 'simulationType':
%   1) 'serial': The ensemble members are simulated serially. This is slow,
%      but can be used for small ensembles/problems and is handy for
%      debugging.
%   2) 'background': The ensemble class runs its members in parallel by
%      spawning new matlab sessions in the background. 
%   3) 'parallel': Parallel execution of the ensemble members by using the
%      Parallel Computing Toolbox. If the toolbox is not found on your
%      system, this option is changed to 'background'.
% The optional parameter 'reset' is used to delete any existing results and
% start the ensemble simulation from scratch.
ensemble = MRSTEnsemble(baseCase, samplesRH, qoiState       , ...
                        'directory'         , dataDir     , ...
                        'simulationStrategy', 'background', ...
                        'reset'             , true       );

%% Simulate the ensemble
% All ensemble members can now be executed with a single line of code and 
% if the 'background' parallelization mode is selected, we show a progress
% bar indicating when the various
close all, ensemble.simulateEnsembleMembers('plotProgress', true);

% It is also possible to run only a subset of the ensemble members by
% calling the function 
%      ensemble.simulationEnsembleMembers(20:30 'plotProress', true);
% for a given range of ensemble members (here ensemble members 20 to 30).


%% Plot results
close all
[x,y] = ndgrid(linspace(0,1000, baseCase.options.ncells));
[s_avg, s_var, s] = ensemble.qoi.getQoIMean();
for i = 1:10:ensembleSize
    baseCase.plot(s{i}.sW); caxis([0,1]);
    title(sprintf('Water saturation, member %d', i));
end
baseCase.plot(s_avg.sW); caxis([0,1]);
title('Water saturation, ensemble avg');

%% References
% [1] "Simulation of Stationary Gaussian Processes in [0, 1]^d", Andrew
% Wood and Grace Chan, Journal of Computational and Graphical Statistics,
% Vol.3 Number 4 (1994), pages 409-432

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
