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

%%
mrstModule add diagnostics sequential linearsolvers
levels = {MLMCLevel(1, 'solver', @(problem) flowDiagnosticsSolver(problem)), ...
          MLMCLevel(2)};
      
%%
dataDir = fullfile(mrstOutputDirectory(), 'ensemble', 'tutorial-mlmc');
mlmc = MLMCSimulator(example, samples, qoi, levels        , ...
                        'directory'         , dataDir     , ...
                        'simulationStrategy', 'serial', ...
                        'verboseSimulation', true);

mlmc.runSimulation('batchSize', 4, 'plotProgress', false, 'minSamples', 10, 'maxSamples', 1000);

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
