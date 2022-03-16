%% Monte carlo and multilevel Monte Carlo simulations
% This tutorial goes through how to set up and run Monte Carlo and
% multilevel Monte Carlo simlations for uncertainty quantification in MRST
% using the using the ensemble module

%% Add necessary modules
mrstModule add ensemble
mrstModule add test-suite
mrstModule add ad-core ad-props ad-blackoil
mrstModule add mrst-gui

%% Set up base case
baseCase = TestCase('qfs_wo', 'ncells', 64);
problem = baseCase.getPackedSimulationProblem();

%% Set up function for generating samples
generatorFn = @(problem, seed) ...
    generateRockSample(baseCase.model.G.cartDims, ...
                                              'seed'    , seed, ...
                                              'toVector', true, ...
                                              'std_perm', 0.5 , ...
                                              'std_poro', 0.1 );
samples = RockSamples('generatorFn', generatorFn);

%% Inspect ensemble subset
% Inspect rock samples
baseCase.plot(problem.SimulatorSetup.model.rock, 'log10', true); colormap(pink);
for i = 1:5
    data    = samples.getSample(i, problem);   % Get sample number i
    problem = samples.setSample(data, problem); % Set sample to problem
    % Inspect rock sample
    baseCase.plot(problem.SimulatorSetup.model.rock, 'log10', true); colormap(pink);
end

%% Define recovery factor
qoi = RecoveryFactorQoI();

%% Set up Monte Carlo simulator
dataDir = fullfile(mrstOutputDirectory(), 'ensemble', 'tutorial-mc', 'mc');
mc = MCSimulator(baseCase, samples, qoi, ...
                        'directory'         , dataDir     , ...
                        'simulationStrategy', 'background', ...
                        'verboseSimulation' , true        );

%% Run monte carlo simulation
batchSize    = 4*maxNumCompThreads(); % Maximum number of samples we 
                                      % simulate in each iteration
relTolerance = 1e-3; % We aim at an RMSE = 1e-3*estimated recovery factor
maxSamples   = batchSize*4; % Maximum number of samples we allow for. This 
                            % must be increased to reach the prescribed
                            % tolerance

% Run simulation. This will bring up three figures: one showing the
% ensemble simulation progress, one showing a histogram of the recovery
% factor estimates, and one illustrating the evolution of the estimate as
% we simulate more samples
close all
mc.runSimulation('batchSize'   , batchSize   , ...
                 'maxSamples'  , maxSamples  , ...
                 'relTolerance', relTolerance);

%% Define coarse levels for multilevel Monte Carlo
mrstModule add coarsegrid

cartDims = baseCase.model.G.cartDims;
n = 3;
partitions = cell(n,1);
for i = 1:n
    partitions{i} = partitionCartGrid(cartDims, cartDims./2^i);
end
partitions = [(1:baseCase.model.G.cells.num)'; partitions];

nc = cumsum([0;cellfun(@max, partitions)]);

partitions0 = partitions;
levels = cell(n+1,1);
levels{1} = MLMCLevel(1);

for i = 2:n+1
    
    partition = partitions0{i};
    
    level = repmat(i, baseCase.model.G.cells.num, 1);
    level([baseCase.schedule.control(1).W.cells]) = 1;
    
    for j = i:-1:2

        cells = accumarray(partition, level < j) > 0;
        cells = cells(partition);

        partition(cells) = partitions0{j-1}(cells);

    end
    
    partition = compressPartition(partition);
    partitions{i} = partition;
    
end

partitions = flipud(partitions);

for i = 1:n
    levels{i} = SpatialUpscalingLevel(i, partitions{i});
end
levels{end} = MLMCLevel(n+1);

%% Set up multilevel Monte Carlo simulator
dataDir = fullfile(mrstOutputDirectory(), 'ensemble', 'tutorial-mc', 'mlmc');
mlmc = MLMCSimulator(baseCase, samples, qoi, levels, ...
                     'directory'         , dataDir     , ...
                     'simulationStrategy', 'background', ...
                     'verboseSimulation' , true        );

%% Show MLMC layers
close all
seed = 1;
cax = [];
for i = mlmc.numLevels():-1:1
    problem = mlmc.levels{i}.levels{end}.getBaseProblem();
    sample  = mlmc.levels{i}.levels{end}.samples.getSample(seed, problem);
    problem = samples.setSample(data, problem); % Set sample to problem
    % Inspect rock sample
    baseCase.figure();
    plotToolbar(problem.SimulatorSetup.model.G, ...
                problem.SimulatorSetup.model.rock, 'log10', true); colormap(pink);
    baseCase.setAxisProperties(gca);
    if isempty(cax)
        ax = gca;
        cmin = min(ax.Children(1).CData);
        cmax = max(ax.Children(1).CData);
        cax = [cmin, cmax];
    end
    caxis(cax);
end

%% Run multilevel Monte Carlo simulation
% Run simulation. This will bring up two figures for each level: one
% showing the ensemble simulation progress, one showing a histogram of the
% level estimate. It will also bring up one figure illustrating the
% evolution of the recovery factor estimate as we simulate more samples
close all
mlmc.runSimulation('batchSize'   , batchSize   , ...
                   'maxSamples'  , maxSamples  , ...
                   'relTolerance', relTolerance);

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
