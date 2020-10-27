classdef MRSTHistoryMatchingEnsemble < MRSTEnsemble
    
    properties
        
        % Most properties inherited from MRSTEnsemble
        
        historyMatchingIteration = 1 
        
        
    end
    
    methods
        %function ensemble = MRSTHistoryMatchingEnsemble(mrstExample, samples, qoi, varargin)
        %    ensemble = ensemble@MRSTEnsemble(mrstExample, samples, qoi, varargin{:});
        %end
        
        %-----------------------------------------------------------------%
        function updatedSample = doHistoryMatching(ensemble)
            
            % Get observations, obs error cov, and scaling factors
            [obs, scaling] = ensemble.qoi.getObservationAndScaling();
            ensembleObs = ensemble.getEnsembleQoI();
            R = ensemble.qoi.getObservationErrorCov();
            
            [obs, ensembleObs, R] = ensemble.applyScaling(obs, ensembleObs, R, scaling);
            [obs, ensembleObs, R] = ensemble.removeObsoleteObservations(obs, ensembleObs, R);
            
            ensembleObs = ensemble.perturbEnsembleQoI(ensembleObs);
            
            ensembleParameters = ensemble.getEnsembleSamples();
            
            analysisParameters = ensemble.enkf(ensembleParameters, obs, ensembleObs, R);
            
            updatedSample = ensemble.samples.setSampleVectors(analysisParameters);            
        end

        
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
        
        function ensembleQoI = perturbEnsembleQoI(ensemble, ensembleQoI)
            % Add unbiased observation error across the ensemle according
            % to the observation error covariance 
            
            obserror = randn(size(ensembleQoI));

            % Avoid biased observation errors by "normalizing" across
            % the ensemble
            obserror = obserror - mean(obserror, 2);
            obserror = obserror ./ std(obserror, 0, 2);

            ensembleQoI = ensembleQoI ...
                + sqrt(ensemble.qoi.getObservationErrorCov())*obserror;
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
        function defaultPath = getDefaultPath(ensemble)
            
            mainFolder = fullfile(mrstOutputDirectory(), 'historyMatching', ensemble.setup.name);
            if ensemble.deleteOldResults && exist(mainFolder, 'dir')
                rmdir(mainFolder, 's');
            end
                        
            defaultPath = fullfile(mrstOutputDirectory(), 'historyMatching', ...
                ensemble.setup.name, num2str(ensemble.historyMatchingIteration));
        end
                
        
    end
end

