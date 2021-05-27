function [estimate, variance, cost] = simulateMC(ensemble, varargin)
    opt = struct('tolerance', 1e-3, ...
                 'batchSize', 10  , ... 
                 'maxSamples' );
    opt = merge_options(opt, varargin{:});
    
    estimate = 0;
    
    while sqrt(variance) > opt.tolerance
        estimate.simulateEnsembleMembers(opt.batchSize)
    end
    
end