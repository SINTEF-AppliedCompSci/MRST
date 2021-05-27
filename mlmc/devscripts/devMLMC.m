mrstModule add ensemble
mrstModule add example-suite
mrstModule add ad-core ad-props ad-blackoil
mrstModule add mrst-gui

%%
example = MRSTExample('qfs_wo');
problem = example.getPackedSimulationProblem();

%%
generatorFn = @(problem, seed) ...
    generateRockSample(problem.SimulatorSetup.model.G.cartDims, ...
                                              'seed'    , seed, ...
                                              'toVector', true, ...
                                              'std_perm', 0.5 , ...
                                              'std_poro', 0.1 );
samples = RockSamples('generatorFn', generatorFn);

%% Inspect ensemble subset
% Inspect rock samples
example.plot(problem.SimulatorSetup.model.rock, 'log10', true); colormap(pink);
for i = 1:5
    data    = samples.getSample(i, problem);   % Get sample number i
    problem = samples.setSample(data, problem); % Set sample to problem
    % Inspect rock sample
    example.plot(problem.SimulatorSetup.model.rock, 'log10', true); colormap(pink);
end

%%
qoi = RecoveryFactorQoI();

%%
dataDir = fullfile(mrstOutputDirectory(), 'ensemble', 'tutorial-mc');
mc = MCSimulator(example, samples, qoi             , ...
                        'directory'         , dataDir     , ...
                        'simulationStrategy', 'background', ...
                        'verboseSimulation', true);

%%
close all
mc.runMonteCarloSimulation('batchSize', 20, 'maxSamples', 10000, 'relTolerance', 1e-4);