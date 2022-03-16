classdef MRSTHistoryMatchingEnsemble < MRSTEnsemble
    
    properties
        
        % Most properties inherited from MRSTEnsemble
        
        historyMatchingIteration = 1
        historyMatchingSubIteration = 1
        esmdaIterations = 1
        
        method = 'EnKF'
        % Allowed options: 
        % 'EnKF': Standard stochastic ensemble Kalman filter/smoother, with
        %         perturbed observations.
        % 'wrongEnKF: Same as 'EnKF', but with observation errors added to
        %         the simulated observations.
        
        mainDirectory; 
        alpha = [1];
        % Folders are organized as follows:
        % mainDirectory/historyMatchingIteration/historyMatchingSubIteration/<ensemble_data>
        
        qoiArchive = {{}};
        
        storeHistoryMatching = true;
        
        originalSchedule
        
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function ensemble = MRSTHistoryMatchingEnsemble(baseCase, samples, qoi, varargin)
            
            opt = struct('prepareSimulation', true, ...
                         'observationProblem', struct([]), ...
                         'perturbObservations', true);
            
            [opt, extra] = merge_options(opt, varargin{:});
            
            ensemble = ensemble@MRSTEnsemble(baseCase, samples, qoi, ...
                                             extra{:}, 'prepareSimulation', false);
            
            % Check that the history-matching-specific input values make
            % sense
            if ensemble.esmdaIterations == 1
                % esmdaIterations not specified, set it according to alpha
                ensemble.esmdaIterations = numel(ensemble.alpha);
            elseif numel(ensemble.alpha) == 1 && ensemble.alpha == 1
                % alpha not specified, set it according to esmdaIterations
                ensemble.alpha = ones(1, ensemble.esmdaIterations)*ensemble.esmdaIterations;
            end
            assert(ensemble.esmdaIterations == numel(ensemble.alpha), ...
                'Number of ES-MDA iterations do not match number of alpha values');
            assert(sum(1./ensemble.alpha) > 0.999 && sum(1./ensemble.alpha) < 1.001, ...
                '1/alpha does not sum to 1');
            
            ensemble.originalSchedule = ensemble.setup.schedule;

            if ensemble.storeHistoryMatching
                if ~exist(ensemble.mainDirectory, 'dir')
                    mkdir(ensemble.mainDirectory);
                end
                save(fullfile(ensemble.mainDirectory, 'ensemble.mat'), 'ensemble');
            end
            
            % if we chose to reset and delete any old results, we need to
            % rebuild the folder structure
            ensemble.getIterationPath();
            
            % Check observations
            if isempty(opt.observationProblem)
               assert(~isempty(ensemble.qoi.observationResultHandler), ...
                   'Missing observations. Provide observations to the QoI directly, or as an observationProblem in the ensemble constructor');
            else
                assert(isempty(ensemble.qoi.observationResultHandler), ...
                    'Observations provided twice. Please only provide observations to the QoI directly, or as an observationProblem in the ensemble constructor. Not both');
                ensemble.qoi = ensemble.qoi.setObservations(opt.observationProblem, ...
                                                            'perturb', opt.perturbObservations);
            end
            
            
            % Prepare ensemble
            if opt.prepareSimulation
                ensemble.prepareEnsembleSimulation();
            end
            
        end
        
        
        function simulateEnsembleMembers(ensemble, varargin)
            opt = struct('plotIntermediateQoI', false);
            [opt, extra] = merge_options(opt, varargin{:});
           
            ensemble.simulateEnsembleMembers@MRSTEnsemble( ...
                'plotIntermediateQoI', opt.plotIntermediateQoI, varargin{:});
            
            if ensemble.storeHistoryMatching
                samples = ensemble.samples;
                save(fullfile(ensemble.directory, 'samples.mat'), 'samples');
            end
            
        end
  
                
        %-----------------------------------------------------------------%
        function doHistoryMatching(ensemble)
            
            for i = 1:ensemble.esmdaIterations
                if i > 1
                    progressTitle = sprintf('Simulating intermediate ensemble %d out of %d', i-1, ensemble.esmdaIterations-1);
                    ensemble.simulateEnsembleMembers('progressTitle', progressTitle);
                end
                
                if ensemble.verbose, fprintf('Starting history matching iteration (%d, %d)\n', ensemble.historyMatchingIteration, ensemble.historyMatchingSubIteration); end
                updatedSamples = ensemble.doHistoryMatchingSingle(ensemble.alpha(i));
                if ensemble.verbose, fprintf('Done history matching iteration (%d, %d)\n', ensemble.historyMatchingIteration, ensemble.historyMatchingSubIteration); end
                
                ensemble.updateHistoryMatchingIterations(updatedSamples);
            end
        end
        
        
        %-----------------------------------------------------------------%
        function updatedSample = doHistoryMatchingSingle(ensemble, alpha)
            
            % Get observations, obs error cov, and scaling factors
            [obs, scaling] = ensemble.qoi.getObservationAndScaling();
            ensembleObs = ensemble.getEnsembleQoI();
            R = ensemble.qoi.getObservationErrorCov();
            
            [obs, ensembleObs, R] = ensemble.applyScaling(obs, ensembleObs, R, scaling);
            [obs, ensembleObs, alphaR] = ensemble.removeObsoleteObservations(obs, ensembleObs, alpha*R);
            
            obsPerturbations = ensemble.perturbEnsembleQoI(ensembleObs, alphaR, 'returnPerturbations', false);
            
            ensembleParameters = ensemble.getEnsembleSamples();
            
            analysisParameters = ensemble.enkf(ensembleParameters, obs, ensembleObs, obsPerturbations, alphaR);
            
            updatedSample = ensemble.samples.setSampleVectors(analysisParameters);            
        end

        %-----------------------------------------------------------------%
        function updateHistoryMatchingIterations(ensemble, samples)
            % This function is called after we have update the ensemble 
            
            % Set new samples to the ensemble class
            ensemble.samples = samples;
            
            % Store QoI
            ensemble.qoiArchive{ensemble.historyMatchingIteration}{ensemble.historyMatchingSubIteration} = ensemble.qoi;
            
            if ensemble.historyMatchingSubIteration == ensemble.esmdaIterations
                ensemble.historyMatchingIteration = ensemble.historyMatchingIteration + 1;
                ensemble.historyMatchingSubIteration = 1;
            else
                ensemble.historyMatchingSubIteration = ensemble.historyMatchingSubIteration + 1;
            end
            
            ensemble.directory = ensemble.getIterationPath();
            
            % Clear the QoI result handler so that it can be recreated
            % pointing to the new result directory
            ensemble.qoi.ResultHandler = []; 
            baseProblem = ensemble.getBaseProblem();
            ensemble.qoi = ensemble.qoi.validateQoI(baseProblem);
            
            ensemble.setupSimulationStatusHandler('');
            
            ensemble.prepareEnsembleSimulation('force', true);
            
        end

        %-----------------------------------------------------------------%
        function xF = enkf(ensemble, xF, obs, obsE, obsPert, R)
           
            % Mixing the notation of the old EnKF module and the review
            % paper by Vetra-Carvalho et al
            Ny = numel(obs);
            Ne = ensemble.num;
    
            S = obsE - mean(obsE, 2);
            
            if strcmp(ensemble.method, 'EnKF')
                obs = repmat(obs, 1, Ne) + obsPert;
            elseif strcmp(ensemble.method, 'wrongEnKF')
                S = S + obsPert;
            end
            
            
            S = S / sqrt(Ne-1);
            
            if Ny <= Ne
                C = S*S' + R;
                W = S'/C;
            else 
                W = (eye(Ne) + (S'/R)*S)  \ S'/R;
            end
            
            xF = xF + (1/sqrt(Ne-1))*(xF - mean(xF, 2))*W*(obs - obsE);
            
        end
        
        
        %-----------------------------------------------------------------%
        function ensembleQoI = getEnsembleQoI(ensemble)
            
            % Check that all ensemble members are computed
            for i = 1:ensemble.num
                assert(ensemble.qoi.isComputed(i), ... 
                    'Ensure that all ensemble members are simulated before calling getEnsembleQoI')
            end
            
            qoivector1 = ensemble.qoi.getQoIVector(1);
            ensembleQoI = zeros(size(qoivector1, 1), ensemble.num);
            ensembleQoI(:,1) = qoivector1;
            
            for i = 2:ensemble.num
                ensembleQoI(:, i) = ensemble.qoi.getQoIVector(i);
            end
        end
        
        %-----------------------------------------------------------------%
        function perturbations = perturbEnsembleQoI(ensemble, ensembleQoI, R, varargin)
            % Sample unbiased observation error across the ensemle according
            % to the observation error covariance. By default this error is
            % added to the ensembleQoI, but can also be returned directly
            
            opt = struct('returnPerturbations', true);
            [opt, extra] = merge_options(opt, varargin{:});
            
            obserror = randn(size(ensembleQoI));

            % Avoid biased observation errors by "normalizing" across
            % the ensemble
            obserror = obserror - mean(obserror, 2);
            obserror = obserror ./ std(obserror, 0, 2);
   
            perturbations = sqrt(R)*obserror;
            if opt.returnPerturbations
                perturbations = ensembleQoI + perturbations;
            end
        end
        
        %-----------------------------------------------------------------%
        function ensembleSamples = getEnsembleSamples(ensemble)
            % Collect all samples in a matrix so that each column represent
            % the sample for a single ensemble member.
            
            ensembleSamples = ensemble.samples.getSampleVectors();
        end
        
        %-----------------------------------------------------------------%
        function [obs, ensembleObs, R] = applyScaling(ensemble, obs, ensembleObs, R, scaling)
            obs = obs./scaling;
            ensembleObs = ensembleObs./scaling;
            R = R./(scaling.^2);
        end
        
        %-----------------------------------------------------------------%
        function [obs, ensembleObs, R] = removeObsoleteObservations(ensemble, obs, ensembleObs, R)
                
            % Following the removal structure as the old EnKF module.
            
            remove1 = find(max(abs(repmat(obs,1,ensemble.num) - ensembleObs),[],2) < eps);
            remove2 = find(std(ensembleObs,[],2) < 10*eps);
            remove3 = find(isnan(obs));
            
            remove  = union(remove1,remove2); 
            remove = union(remove,remove3);
            
            indicesToKeep = setdiff(1:length(obs),remove);
            
            obs = obs(indicesToKeep);
            ensembleObs = ensembleObs(indicesToKeep, :);
            R = R(indicesToKeep, indicesToKeep);
        end
        
        %-----------------------------------------------------------------%
        function updateHistoryMatchingInterval(ensemble, historyMatchDtRange)
            % Takes a range of step indices as input
            ensemble.qoi.historyMatchDtRange = historyMatchDtRange;
            ensemble.qoi.dt = ensemble.originalSchedule.step.val(1:historyMatchDtRange(end));
            
            ensemble.setup.schedule.step.val = ensemble.originalSchedule.step.val(1:historyMatchDtRange(end));
            ensemble.setup.schedule.step.control = ensemble.originalSchedule.step.control(1:historyMatchDtRange(end));

            ensemble.prepareEnsembleSimulation('force', true);
        end
        
        %-----------------------------------------------------------------%
        function defaultPath = getDefaultPath(ensemble)
            
            defaultPath = fullfile(mrstOutputDirectory(), 'historyMatching', ensemble.setup.name);
            %if ensemble.reset && exist(ensemble.mainDirectory, 'dir')
            %    rmdir(ensemble.mainDirectory, 's');
            %end
        end
        
        %-----------------------------------------------------------------%
        function iterationDataPath = getIterationPath(ensemble)
            mainIterationDataPath = fullfile(ensemble.mainDirectory, ... 
                                             num2str(ensemble.historyMatchingIteration));
            iterationDataPath = fullfile(mainIterationDataPath, ...
                                         num2str(ensemble.historyMatchingSubIteration));
            
            if ~exist(mainIterationDataPath, 'dir')
                mkdir(mainIterationDataPath);
            end
            if ~exist(iterationDataPath, 'dir')
                mkdir(iterationDataPath);
            end
        end
                
        %-----------------------------------------------------------------%
        function reset(ensemble, varargin)
            % Reset the simulation by simply removing the entire result
            % folder (and subfolders) for the case (including iterations
            % and sub-iterations).
            
            % TODO: Align it with reset@MRSTEnsemble somehow...
            % TODO: Do reset based on opt.prepareSimulation...
            
            opt = struct('prompt', true, ...
                         'prepareSimulation', false);
            opt = merge_options(opt, varargin{:});
            if ~exist(ensemble.mainDirectory, 'dir'), return, end
            if opt.prompt
                prompt = sprintf(['Delete all data for %s? y/n [n]: '], ...
                                              ensemble.setup.name     );
                if ~strcmpi(input(prompt, 's'), 'y')
                    fprintf('Ok, will not remove files.\n');
                    return
                end
            end
            rmdir(ensemble.mainDirectory, 's');
        end     
    
        %-----------------------------------------------------------------%
        function h = plotQoI(ensemble, varargin)
            % Creates plot(s) of the quantity of interest for all simulated
            % ensemble members
            %
            % SYNOPSIS:
            %   h = ensemble.plotQoI();
            %
            % OPTIONAL PARAMETERS:
            %   'h' - Figure handle 
            opt = struct('h',                   [], ...
                         'color',               [], ... % color input will for now be ignored
                         'subIterations',       false, ... 
                         'iterations',          true,  ...
                         'plotObservation',     true,  ...
                         'plotTruth',           false, ...
                         'cmapName',            'tatarizeMap', ...
                         'cmapRealizationsMin', 6, ...
                         'legend',              {{}}, ...
                         'savefig',             false, ...
                         'savefolder',          '');
            [opt, extra] = merge_options(opt, varargin{:});
            
            ci = 1;
            
            numPlots = 1 + opt.iterations*(numel(ensemble.qoiArchive) + ...
                                           opt.subIterations*(ensemble.esmdaIterations - 1));
            cmap = feval(opt.cmapName, max(numPlots, opt.cmapRealizationsMin));
            
            alreadyOpenFigures = length(findobj('type', 'figure'));
            
            h = opt.h;
            if opt.iterations && ~isempty(ensemble.qoiArchive{1})
                for i = 1:numel(ensemble.qoiArchive)
                    for j = 1:ensemble.esmdaIterations 
                        if opt.subIterations || j == 1 
                            h = ensemble.qoiArchive{i}{j}.plotEnsembleQoI(ensemble, h, ...
                                                                          'color', cmap(ci, :), ...
                                                                          'plotObservation', false, ...
                                                                          'plotTruth', false, ...
                                                                          'alreadyOpenFigures', alreadyOpenFigures, ...
                                                                          'savefig', false, ...
                                                                          extra{:});
                            ci = ci + 1;
                        end
                    end
                end
            end
            h = ensemble.qoi.plotEnsembleQoI(ensemble, h, ...
                                             'color', cmap(ci, :), ...
                                             'plotObservation', opt.plotObservation, ...
                                             'plotTruth', opt.plotTruth, ...
                                             'legend', opt.legend, ...
                                             'alreadyOpenFigures', alreadyOpenFigures, ...
                                             'savefig', opt.savefig, ...
                                             extra{:});                               


        
        end
        
        
    end % public methods
    
    methods (Access = protected)
        
        %-----------------------------------------------------------------%
       function setUpDirectory(ensemble)
            % Set up directory name correctly.

            if isempty(ensemble.directory)
                ensemble.directory = ensemble.getDefaultPath();
            end
            ensemble.mainDirectory = ensemble.directory;
            ensemble.directory = ensemble.getIterationPath();
       end
       
    end
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
