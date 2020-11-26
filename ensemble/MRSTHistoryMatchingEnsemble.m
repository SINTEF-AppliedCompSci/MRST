classdef MRSTHistoryMatchingEnsemble < MRSTEnsemble
    
    properties
        
        % Most properties inherited from MRSTEnsemble
        
        historyMatchingIteration = 1
        historyMatchingSubIteration = 1
        esmdaIterations = 1
        
        mainDirectory; 
        alpha = [1];
        % Folders are organized as follows:
        % mainDirectory/historyMatchingIteration/historyMatchingSubIteration/<ensemble_data>
        
        qoiArchive = {{}};
        
        originalSchedule
        
    end
    
    methods
        
        %-----------------------------------------------------------------%
        %function ensemble = MRSTHistoryMatchingEnsemble(mrstExample, samples, qoi, varargin)
        % 
        %    ensemble = ensemble@MRSTEnsemble(mrstExample, samples, qoi, varargin{:});
        
        function midConstructor(ensemble)
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
        end
        
  
                
        %-----------------------------------------------------------------%
        function doHistoryMatching(ensemble)
            
            for i = 1:ensemble.esmdaIterations
                if i > 1
                    ensemble.simulateEnsembleMembers();
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
            
            ensembleObs = ensemble.perturbEnsembleQoI(ensembleObs, alphaR);
            
            ensembleParameters = ensemble.getEnsembleSamples();
            
            analysisParameters = ensemble.enkf(ensembleParameters, obs, ensembleObs, alphaR);
            
            updatedSample = ensemble.samples.setSampleVectors(analysisParameters);            
        end

        %-----------------------------------------------------------------%
        function updateHistoryMatchingIterations(ensemble, updatedSamples)
            % This function is called after we have update the ensemble 
            
            ensemble.samples = updatedSamples;
            
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
            
            ensemble.prepareEnsembleSimulation('force', true);
            
            %if numel(ensemble.spmdEnsemble) > 0 && strcmp(ensemble.simulationStrategy, 'parallel')
            %    % Call recursively on spmdEnsembles (if any)
            %    spmdEns = ensemble.spmdEnsemble;
            %    spmd
            %        spmdEns.updateHistoryMatchingIterations(updatedSamples);
            %    end
            %end
        end

        %-----------------------------------------------------------------%
        function xF = enkf(ensemble, xF, obs, obsE, R)
           
            % Mixing the notation of the old EnKF module and the review
            % paper by Vetra-Carvalho et al

            S = obsE - mean(obsE, 2);
            
            Ny = numel(obs);
            Ne = ensemble.num;
            if Ny <= Ne
                C = S*S'/(Ne-1) + R;
                W = S'/C;
            else 
                W = (eye(Ne) + (S'/R)*S)  \ S'/R;
            end
            
            xF = xF + (1/(Ne-1))*(xF - mean(xF, 2))*W*(obs - obsE);
            
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
        function ensembleQoI = perturbEnsembleQoI(ensemble, ensembleQoI, R)
            % Add unbiased observation error across the ensemle according
            % to the observation error covariance 
            
            obserror = randn(size(ensembleQoI));

            % Avoid biased observation errors by "normalizing" across
            % the ensemble
            obserror = obserror - mean(obserror, 2);
            obserror = obserror ./ std(obserror, 0, 2);

            ensembleQoI = ensembleQoI ...
                + sqrt(R)*obserror;
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
            remove2 = find(std(ensembleObs,[],2) < eps);
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
            % takes a range of step 
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
            iterationDataPath = fullfile(ensemble.mainDirectory, ... 
                                         num2str(ensemble.historyMatchingIteration), ...
                                         num2str(ensemble.historyMatchingSubIteration));
        end
                
        %-----------------------------------------------------------------%
        function reset(ensemble, varargin)
            % Reset the simulation by simply removing the entire result
            % folder (and subfolders) for the case (including iterations
            % and sub-iterations).
            
            % TODO: Align it with reset@MRSTEnsemble somehow...
            
            opt = struct('prompt', true);
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
                         'legend',              {{}});
            [opt, extra] = merge_options(opt, varargin{:});
            
            ci = 1;
            
            numPlots = 1 + opt.iterations*(numel(ensemble.qoiArchive) + ...
                                           opt.subIterations*(ensemble.esmdaIterations - 1));
            cmap = feval(opt.cmapName, max(numPlots, opt.cmapRealizationsMin));
            
            h = opt.h;
            if opt.iterations
                for i = 1:numel(ensemble.qoiArchive)
                    for j = 1:ensemble.esmdaIterations 
                        if opt.subIterations || j == ensemble.esmdaIterations 
                            h = ensemble.qoiArchive{i}{j}.plotEnsembleQoI(ensemble, h, ...
                                                                          'color', cmap(ci, :), ...
                                                                          'plotObservation', false, ...
                                                                          'plotTruth', false, ...
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


