%% History matching unknown but constant permeability
% in a square 2D domain with a single producer/injector pair in opposing
% corners, both controlled by rates with absolute value 1.
% 
mrstModule add ad-core ad-blackoil mrst-gui ad-props ...
    example-suite incomp ensemble 

mrstVerbose off

tic

%% Some folder config

useKozenyCarman = false;

topDirectory = fullfile(mrstOutputDirectory(), 'PINN_examples', 'twoProds_rockSample_freeRock');
if useKozenyCarman
    topDirectory = fullfile(mrstOutputDirectory(), 'PINN_examples', 'twoProds_rockSample_KC');
end
    
referenceDirectory = fullfile(topDirectory, 'truth');
historyMatchingDirectory = fullfile(topDirectory, 'history_matching');

%% Setup full reference model
rng(0);
referenceExample = MRSTExample('pinn_da_case_two_producers', ...
                               'nsteps', 100, ...
                               'tend', 3.0);


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

obsStdDevFlux = 0.01;
obsStdDevBhp = 0.01;
obsStdDev = [obsStdDevFlux, obsStdDevFlux, obsStdDevBhp];

qoi = WellQoIHM('wellNames', {referenceExample.schedule.control.W.name}, ...
                'names', {'qOs', 'qWs', 'bhp'}, ...
                'observationCov', obsStdDev.^2);
                %'observationResultHandler', observationQoI.ResultHandler, ...
                %'truthResultHandler', truthResultHandler, ...
%% Samples

ensembleSize = 100;

rng(100);
rockData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    rockData{i}.poro = 1; %repmat(1, referenceExample.options.ncells^2, 1);
    %rockData{i}.perm = randn(1)*0.05 + 1.1; %repmat(randn(1)*0.05 + 1, referenceExample.options.ncells^2, 1);
    rockData{i}.perm = exp(2*randn(1,1));

    % The reference model is:
    % p = gaussianField(G.cartDims, [0.2 0.4], [11 3], 2.5);

    p = gaussianField(referenceExample.model.G.cartDims, [0.1 0.6], [11 3], 7.5);
    K = p.^3.*(1e-5)^2./(0.81*72*(1-p).^2);
    
    rockData{i}.poro = p;
    rockData{i}.perm = K;
end

rockSamples = RockSamplesHM('data', rockData, ...
                            'minPoroValue', 0.05, ...
                            'maxPoroValue', 0.95, ...
                            'minPermValue', 5*eps, ...
                            'maxPermValue', 1e-8, ...
                            'useKozenyCarman', useKozenyCarman);
% useKozenyCarman = false works good!


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
observationIndices = (5:5:totalNumberOfTimesteps);

%%
ensemble.simulateEnsembleMembers('plotIntermediateQoI', false);

%%
if true
    hm_legend = {'observations', 'truth', 'prior mean'};

    close all
    disp('Plotting prior in progress, please wait...');
    ensemble.plotQoI('subplots', false, 'clearFigure', false, ...
        'cmapName', 'lines', ...
        'plotTruth', true, ...
        'subIterations', false, ...
        'observationIndices', observationIndices, ...
        'legend', hm_legend); %, ...
    %    'savefig', true, ...
    %    'saveFolder', fullfile(historyMatchingDirectory, 'prior_figures'));
    disp('Plotting prior completed');
end

%%

ensemble.updateHistoryMatchingInterval(observationIndices);

ensemble.doHistoryMatching()

%%
ensemble.updateHistoryMatchingInterval(1:totalNumberOfTimesteps);

ensemble.simulateEnsembleMembers('plotIntermediateQoI', false);



%%

plotIntermediate = false;
hm_legend = {'observations', 'truth', 'posterior mean', 'prior mean'};
if plotIntermediate
    hm_legend = {'observations', 'truth', 'posterior mean', 'ES-MDA it 3',...
                   'ES-MDA it 2', 'ES-MDA it 1', 'prior mean'}; %#ok
end
    
    
close all
disp('Plotting posterior in progress, please wait...');
ensemble.plotQoI('subplots', false, 'clearFigure', false, ...
    'cmapName', 'lines', ...
    'plotTruth', true, ...
    'subIterations', plotIntermediate, ...
    'observationIndices', observationIndices, ...
    'legend', hm_legend, ...
    'savefig', true, ...
    'saveFolder', fullfile(historyMatchingDirectory, 'posterior_figures'))
disp('Plotting posterior completed');



%%
meanSample = ensemble.samples.getMeanSample();
varSample = ensemble.samples.getVarianceSample();
G = ensemble.setup.model.G;

refPoro = referenceExample.model.rock.poro(:);
refPerm = log(referenceExample.model.rock.perm(:));
postPoro = meanSample.data{1}.poro(:);
postPerm = log(meanSample.data{1}.perm(:));

varPoro = varSample.data{1}.poro(:);
varPerm = varSample.data{1}.perm(:);

minPerm = min(min(refPerm(:)), min(postPerm));
minPoro = min(min(refPoro(:)), min(postPoro));
maxPerm = max(max(refPerm(:)), max(postPerm))-2;
maxPoro = max(max(refPoro(:)), max(postPoro));

figure()
subplot(3,2,1)
plotCellData(G, refPoro)
caxis([minPoro, maxPoro])
title('reference poro')
colorbar();

subplot(3,2,2)
plotCellData(G, refPerm)
caxis([minPerm, maxPerm])
title('reference log(perm)')
colorbar();

subplot(3,2,3)
plotCellData(G, postPoro)
caxis([minPoro, maxPoro])
title('posterior mean poro')
colorbar();

subplot(3,2,4)
plotCellData(G, postPerm)
caxis([minPerm, maxPerm])
title('posterior log(mean perm)')
colorbar();

subplot(3,2,5)
plotCellData(G, sqrt(varPoro))
title('posterior std.dev. poro')
colorbar();

subplot(3,2,6)
plotCellData(G, log(sqrt(varPerm)))
title('posterior log (std.dev. perm)')
colorbar();

figHandle = gcf;
savefilename = fullfile(historyMatchingDirectory, 'parameter_maps.png');
saveas(figHandle, savefilename, 'png');

%% save configuration
if true
    % topDirectory
    %csvwrite(fullfile(topDirectory, 'reference_perm.csv'),  referenceExample.model.rock.perm(:));
    %csvwrite(fullfile(topDirectory, 'reference_poro.csv'),  referenceExample.model.rock.poro(:));
    
    writematrix(referenceExample.model.rock.perm(:), fullfile(topDirectory, 'reference_perm.csv'));
    writematrix(referenceExample.model.rock.poro(:), fullfile(topDirectory, 'reference_poro.csv'));
    
    observationMatrix = [
        ensemble.qoi.observationResultHandler{1}(1).qWs ...
        ensemble.qoi.observationResultHandler{1}(2).qWs, ...
        ensemble.qoi.observationResultHandler{1}(3).bhp
        ];
    writematrix(observationMatrix, fullfile(topDirectory, 'observations.csv'));
    
    truthMatrix = [
        ensemble.qoi.truthResultHandler{1}(1).qWs ...
        ensemble.qoi.truthResultHandler{1}(2).qWs, ...
        ensemble.qoi.truthResultHandler{1}(3).bhp
        ];
    writematrix(truthMatrix, fullfile(topDirectory, 'trueValues.csv'));
    
    
    
end


%%
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
