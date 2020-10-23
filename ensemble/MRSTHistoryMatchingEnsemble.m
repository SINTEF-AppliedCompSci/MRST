classdef MRSTHistoryMatchingEnsemble < MRSTEnsemble
    
    properties
        % Inherited from MRSTEnsemble only
    end
    
    methods
        function MRSTHistoryMatchingEnsemble()
        
        
        
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
        
        %function simulated
    end
end

