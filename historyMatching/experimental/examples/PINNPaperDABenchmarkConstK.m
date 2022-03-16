%% History matching unknown but constant permeability
% in a square 2D domain with a single producer/injector pair in opposing
% corners, both controlled by rates with absolute value 1.
% 
mrstModule add ad-core ad-blackoil mrst-gui ad-props ...
    example-suite incomp ensemble 

mrstVerbose off


%% Some folder config
topDirectory = fullfile(mrstOutputDirectory(), 'PINN_examples', 'const_perm');

referenceDirectory = fullfile(topDirectory, 'truth');
historyMatchingDirectory = fullfile(topDirectory, 'history_matching');

%% Setup full reference model
referenceExample = MRSTExample('pinn_da_case_const_perm', ...
                               'nsteps', 100, ...
                               'tend', 10.0);


%% Pack and run reference problem
referenceProblem = referenceExample.getPackedSimulationProblem('Directory', referenceDirectory);


rerunReferenceModel = true;
plotReferenceModel = true;

if rerunReferenceModel
    clearPackedSimulatorOutput(referenceProblem, 'prompt',  false);
end
simulatePackedProblem(referenceProblem);

if plotReferenceModel
    [refWellSols, refStates, refReports] = getPackedSimulatorOutput(referenceProblem);
    referenceExample.plot(refStates);
    plotWellSols(refWellSols);
end




            
%% Reservoir state qoi 


% We are interested in the pressure in the well cells 
%qoi = ReservoirStateQoI('names', 'pressure', ...
%                        'cells', [1 referenceExample.model.G.cells.num]);

%qoi = WellQoI(...
%    'wellNames', {'P1', 'I1'}, ...
%    'names', {'bhp'}, ...
%    'cumulative', false);

%obsStdDev = 0.01;
obsStdDev = 1;


qoi = WellQoIHM('wellNames', {'I1'}, ...
                'names', {'bhp'}, ...
                'observationCov', obsStdDev^2);
%qoi = WellQoIHM('wellNames', {'P1', 'I1'}, ...
%                'names', {'bhp', 'qWs'}, ...
%                'observationCov', obsStdDev^2);
                %'observationResultHandler', observationQoI.ResultHandler, ...
                %'truthResultHandler', truthResultHandler, ...
%% Samples

ensembleSize = 100;

rockData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    rockData{i}.poro = 1; %repmat(1, referenceExample.options.ncells^2, 1);
    %rockData{i}.perm = randn(1)*0.05 + 1.1; %repmat(randn(1)*0.05 + 1, referenceExample.options.ncells^2, 1);
    rockData{i}.perm = exp(1 + 2*randn(1,1));
end

rockSamples = RockSamplesHM('data', rockData);



%% Create ensemble and run

alpha = [28/3, 7, 4, 2];
%alpha = [1];

ensemble = MRSTHistoryMatchingEnsemble(referenceExample, rockSamples, qoi, ...
    'observationProblem', referenceProblem, ...
    'alpha', alpha, ...
    'directory', historyMatchingDirectory, ...
    'simulationStrategy', 'spmd', ...
    'maxWorkers', 8, ...
    'reset', true, ...
    'verbose', true, ...
    'verboseSimulation', false);

%ensemble.qoi = ensemble.qoi.setObservations(referenceProblem);


totalNumberOfTimesteps = numel(ensemble.originalSchedule.step.val);
%observationIndices = (2:2:floor(totalNumberOfTimesteps/2));
observationIndices = (1:5);
%observationIndices = 1;

%%
ensemble.simulateEnsembleMembers('plotIntermediateQoI', false);

% hm_legend = {'observations', 'truth', 'prior mean'};
% 
% close all
% disp('Plotting prior in progress, please wait...');
% ensemble.plotQoI('subplots', false, 'clearFigure', false, ...
%     'cmapName', 'lines', ...
%     'plotTruth', true, ...
%     'subIterations', false, ...
%     'observationIndices', observationIndices, ...
%     'legend', hm_legend, ...
%     'savefig', true, ...
%     'saveFolder', fullfile(historyMatchingDirectory, 'prior_figures'));
% disp('Plotting prior completed');


%%

ensemble.updateHistoryMatchingInterval(observationIndices);

ensemble.doHistoryMatching()

%%
ensemble.updateHistoryMatchingInterval(1:totalNumberOfTimesteps);

ensemble.simulateEnsembleMembers('plotIntermediateQoI', false);



%%

hm_legend = {'observations', 'truth', 'posterior mean', 'prior mean'};

close all
disp('Plotting posterior in progress, please wait...');
ensemble.plotQoI('subplots', false, 'clearFigure', false, ...
    'cmapName', 'lines', ...
    'plotTruth', true, ...
    'subIterations', false, ...
    'observationIndices', observationIndices, ...
    'legend', hm_legend); %, ...
    %'savefig', true, ...
    %'saveFolder', fullfile(historyMatchingDirectory, 'posterior_figures'))
disp('Plotting posterior completed');



%%

%ensemble.sample
priorSample     = load(fullfile(ensemble.mainDirectory, '1', '1', 'samples.mat'));
posteriorSample = load(fullfile(ensemble.mainDirectory, '2', '1', 'samples.mat'));

% Extract perms
priorSample.samples.transformSampleVectors = false;
posteriorSample.samples.transformSampleVectors = false;

priorPerm = priorSample.samples.getSampleVectors();
posteriorPerm = posteriorSample.samples.getSampleVectors();

priorPerm = priorPerm(1,:);
posteriorPerm = posteriorPerm(1,:);

figure()
histbuckets = 40;
histogram(log(priorPerm), histbuckets);
hold on
histogram(log(posteriorPerm), histbuckets);
plot(log([1 1]), [0 25]);
legend('prior', 'posterior', 'truth')

title('Distribution of log(K)');




%%
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
