classdef MLMCSimulator < MCSimulator
    
    properties
        levels
    end
   
    methods
        %-----------------------------------------------------------------%
        function mlmc = MLMCSimulator(setup, samples, qoi, levels, varargin)
            mlmc = mlmc@MCSimulator(setup, samples, qoi, varargin{:});
            mlmc.levels = cell(numel(levels), 1);
            for i = 1:mlmc.numLevels
                lix = max(i-1,1):i;
                mlmc.levels{i} = MCLevelSimulator(setup, samples, qoi, levels(lix), ...
                                                                      varargin{:});
            end
        end
        
        %-----------------------------------------------------------------%
        function range = runBatch(mlmc, varargin)
            opt = struct('batchSize', []);
            [opt, extra] = merge_options(opt, varargin{:});
            range = nan(sum(opt.batchSize), 1);
            for i = 1:mlmc.numLevels
                if opt.batchSize(i) == 0, continue; end
                maxId = mlmc.getLargestSeed();
                lrange = (1:opt.batchSize(i)) + maxId;
                mlmc.levels{i}.runBatch('range', lrange, extra{:});
                range((1:opt.batchSize(i)) + sum(opt.batchSize(1:i-1))) = lrange;
            end
        end
        
        %-----------------------------------------------------------------%
        function updateStatistics(mlmc, range)
            for i = 1:mlmc.numLevels()
                ids = mlmc.levels{i}.qoi.ResultHandler.getValidIds();
                lrange = range(ismember(range, ids));
                mlmc.levels{i}.updateStatistics(lrange);
            end
            % Get level statistices
            stat = mlmc.getLevelStatistics();
            % Compute statistics
            mlmc.estimate   = sum(stat.estimate);
            mlmc.variance   = sum(stat.variance./stat.numSamples);
            mlmc.rmse       = sqrt(mlmc.variance);
            mlmc.cost       = sum(stat.numSamples.*stat.cost);
            mlmc.numSamples = sum(stat.numSamples);
            % Update history
            mlmc.history{end+1} = struct('estimate'  , mlmc.estimate  , ...
                                         'variance'  , mlmc.variance  , ...
                                         'cost'      , mlmc.cost      , ...
                                         'numSamples', mlmc.numSamples, ...
                                         'rmse'      , mlmc.rmse      );
        end
        
        %-----------------------------------------------------------------%
        function maxId = getLargestSeed(mc)
            maxId = cellfun(@(level) max(level.qoi.ResultHandler.getValidIds()), ...
                                             mc.levels, 'UniformOutput', false);
            maxId = cell2mat(maxId(~cellfun(@isempty, maxId)));
            maxId = max(maxId);
            if isempty(maxId), maxId = 0; end
        end
        
        %-----------------------------------------------------------------%
        function range = getComputedSampleRange(mlmc)
            range = [];
            for i = 1:mlmc.numLevels
                range = [range, ...
                       mlmc.levels{i}.qoi.ResultHandler.getValidIds()]; %#ok
            end
        end
        
        %-----------------------------------------------------------------%
        function reset(mlmc, varargin)
            for i = 1:mlmc.numLevels
                mlmc.levels{i}.reset(varargin{:});
            end
        end
        
        %-----------------------------------------------------------------%
        function n = computeNumSamples(mlmc, tolerance, opt)
            stat = mlmc.getLevelStatistics();
            n = tolerance.^(-2).*sum(sqrt(stat.variance.*stat.cost)) ...
                                   .*sqrt(stat.variance./stat.cost);
            n = ceil(n);
            n = min(n, opt.maxSamples - mlmc.numSamples);
            n = max(n, opt.minSamples - mlmc.numSamples);
            n = min(n, opt.batchSize);
            n = max(n, 0);
        end
        
        %-----------------------------------------------------------------%
        function n = numLevels(mlmc)
            n = numel(mlmc.levels);
        end
        
        function stat = getLevelStatistics(mlmc)
            lstat = cellfun(@(mcl) mcl.getStatistics, mlmc.levels, 'UniformOutput', false);
            stat = struct();
            for fn = {'estimate', 'variance', 'cost', 'numSamples'}
                stat.(fn{1}) = cellfun(@(lstat) lstat.(fn{1}), lstat);
            end
        end
    end
   
end