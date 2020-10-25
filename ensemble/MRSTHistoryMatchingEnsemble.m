classdef MRSTHistoryMatchingEnsemble < MRSTEnsemble
    
    properties
        
        % Most properties inherited from MRSTEnsemble
        
        historyMatchingIteration = 1 
        
        
    end
    
    methods
        %function ensemble = MRSTHistoryMatchingEnsemble(mrstExample, samples, qoi, varargin)
        %    ensemble = ensemble@MRSTEnsemble(mrstExample, samples, qoi, varargin{:});
        %end
        
        function ensembleQoI = getEnsembleQoI(ensemble, varargin)
            opt = struct('perturb', false);
            [opt, extra] = merge_options(opt, varargin{:});
            
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
            
            % Add unbiased observations
            if opt.perturb
                obserror = randn(size(ensembleQoI));
                
                % Avoid biased observation errors by "normalizing" across
                % the ensemble
                obserror = obserror - mean(obserror, 2);
                obserror = obserror ./ std(obserror, 0, 2);
                
                ensembleQoI = ensembleQoI ...
                    + sqrt(ensemble.qoi.getObservationErrorCov())*obserror;
            end
        end
        
        function ensembleSamples = getEnsembleSamples(ensemble)
            % Collect all samples in a matrix so that each column represent
            % the sample for a single ensemble member.
            
            ensembleSamples = ensemble.samples.getSampleVectors();
        end
        
        %-----------------------------------------------------------------%
        function defaultPath = getDefaultPath(ensemble)
            defaultPath = fullfile(mrstOutputDirectory(), 'historyMatching', ...
                ensemble.setup.name, num2str(ensemble.historyMatchingIteration));
        end
                
        
    end
end

