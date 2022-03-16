%% History mathing of 1D reservoir simulation
% In this example, we use well production data to estimate 
mrstModule add ad-core ad-blackoil mrst-gui ad-props ...
    test-suite incomp ensemble

mrstVerbose off


%% Set up and simulate the true solution
% We will here use an identical twin experiment, where we use the same
% problem for both generating the truth and as a base for our ensemble.

rng(0);

trueProblemName = 'ensemble_base_problem_1d_reservoir';
numCells = 10;
trueProblemOptions = {'ncells', numCells, ...
                      'rngseed', 1};
                  
directoryTruth = fullfile(mrstOutputDirectory(), ...
                          'historyMatching', 'truth', ...
                          trueProblemName);
                  
trueCase = TestCase(trueProblemName, trueProblemOptions{:});
trueProblem = trueCase.getPackedSimulationProblem('Directory', directoryTruth);

plotExample = true;
rerunTrueProblemFromScratch = true;
overwriteObservation = rerunTrueProblemFromScratch || false;

if rerunTrueProblemFromScratch
    clearPackedSimulatorOutput(trueProblem, 'prompt', false);
end
simulatePackedProblem(trueProblem);
if plotExample
    [wellSols, states, reports] = getPackedSimulatorOutput(trueProblem);
    trueCase.plot();
    trueCase.plot(states);
end


%% Select and populate samples for the stochastic components in the ensemble

ensembleSize = 40;

configData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    configData{i}.poro = gaussianField(trueCase.model.G.cartDims, [0.2 0.6]); 
    configData{i}.perm = configData{i}.poro.^3.*(1e-5)^2./(0.81*72*(1-configData{i}.poro).^2);
end

samples = RockSamplesHM('data', configData, ...
                        'minPoroValue', 0.02)

%% Select quantity of interest class matching the what we have as observations
% We validate the QoI with the trueProblem, since this will be our ensemble
% base problem as well.
obsStdDev = 0.00005;

qoi = WellQoIHM('wellNames', {'P1'}, 'names', {'qOs'}, ...
                  ... % 'observationResultHandler', observationResultHandler, ...
                  ... % 'truthResultHandler', truthResultHandler, ...
                  'observationCov', obsStdDev^2)

%% Create the ensemble

ensemble = MRSTHistoryMatchingEnsemble(trueCase, samples, qoi, ... 
    'observationProblem', trueProblem, ...
    ... %'directory', uniqueDirectory, ...
    'alpha', [28/3 7 4 2], ...
    'simulationStrategy', 'background', ...
    'maxWorkers', 4, ...
    'verbose', true, ...
    'reset', true...
    )


%% Displaying the observations and observation error cov through the ensemble
disp('observation and scaling vector')
[obsVector, obsScaling] = ensemble.qoi.getObservationAndScaling();
disp('observation error covariance matrix')
ensemble.qoi.getObservationErrorCov()


%% Run ensemble
ensemble.simulateEnsembleMembers();

%% Plot prior
ensemble.plotQoI('subplots', true, 'clearFigure', false, ...
    'cmapName', 'lines', ...
    'plotTruth', true, ...
    'legend', {'observations', 'truth', 'posterior mean'});

%%
perm_range = [1:10];
poro_range = [11:20];
param_ranges = {perm_range, poro_range};

priorDir = fullfile(ensemble.mainDirectory, 'prior');
mkdir(priorDir);
samplesStructure = postProcessRockSamples(ensemble.mainDirectory, param_ranges, ...
    'saveFigFolder', priorDir);

%% Get simulated observations
disp('simulated observations')
size(ensemble.getEnsembleQoI())

%% Get the matrix of ensemble samples 
size(ensemble.getEnsembleSamples())

%% Do history matching
disp('updated sample object:')
ensemble.doHistoryMatching()


%% Run new ensemble
ensemble.simulateEnsembleMembers();

%% Plot 
if numel(ensemble.alpha) == 1
    ensemble.plotQoI('subplots', true, 'clearFigure', false, ...
        'cmapName', 'lines', ...
        'plotTruth', true, ...
        'legend', {'observations', 'truth', 'posterior mean', 'prior mean'});
else
    ensemble.plotQoI('subplots', true, 'clearFigure', false, ...
        'cmapName', 'lines', ...
        'plotTruth', true, ...
        'subIterations', true, ...
        'legend', {'observations', 'truth', 'posterior mean', 'ES-MDA it 3',...
                   'ES-MDA it 2', 'ES-MDA it 1', 'prior mean'});
end
           


%% Plot diff between the qois
%figure
%hold on
%for i = 1:ensemble.num
%    plot(ensemble.qoi.ResultHandler{i}{1}{1} - ensemble.qoiArchive{1}{1}.ResultHandler{i}{1}{1})
%end
%title('difference in production rates')

%% Plot parameter distribution with trueValues
posteriorDir = fullfile(ensemble.mainDirectory, 'posterior');
mkdir(posteriorDir);

samplesStructure = postProcessRockSamples(ensemble.mainDirectory, param_ranges, ...
    'trueRock', trueProblem.SimulatorSetup.model.rock, ...
    'saveFigFolder', posteriorDir);




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
