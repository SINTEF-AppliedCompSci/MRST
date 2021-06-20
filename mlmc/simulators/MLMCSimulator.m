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
        function rangeStat = updateStatistics(mlmc, range)
            rangeStat = cell(mlmc.numLevels,1);
            for i = 1:mlmc.numLevels()
                ids = mlmc.levels{i}.qoi.ResultHandler.getValidIds();
                lrange = range(ismember(range, ids));
                if ~isempty(lrange)
                    rangeStat{i} = mlmc.levels{i}.updateStatistics(lrange);
                end
            end
            % Get level statistices
            stat = mlmc.getLevelStatistics();
            % Compute statistics
            mlmc.numSamples = sum(stat.numSamples);
            mlmc.estimate   = sum(stat.estimate);
            mlmc.variance   = sum(stat.variance./stat.numSamples);
            mlmc.rmse       = sqrt(mlmc.variance);
            mlmc.cost       = sum(stat.numSamples.*stat.cost)./mlmc.numSamples;
            % Update history
            mlmc.history{end+1} = struct('numSamples', mlmc.numSamples, ...
                                         'estimate'  , mlmc.estimate  , ...
                                         'variance'  , mlmc.variance  , ...
                                         'cost'      , mlmc.cost      , ...
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
        function flag = reset(mlmc, varargin)
            flag = true;
            opt = struct('prompt', true);
            [opt, extra] = merge_options(opt, varargin{:});
            dataPath = mlmc.getDataPath();
            if ~exist(dataPath, 'dir'), return, end
            if opt.prompt
                prompt = sprintf(['Delete all data for all levels of' , ...
                                  ' MLMC simulation %s? (sample data' , ...
                                  ' will not be deleted) y/n [n]: '  ], ...
                                                             mlmc.name);
                if ~strcmpi(input(prompt, 's'), 'y')
                    fprintf('Ok, will not remove files.\n');
                    flag = false;
                    return
                end
            end
            for i = 1:mlmc.numLevels
                mlmc.levels{i}.reset('prompt', false, extra{:});
            end
        end
        
        %-----------------------------------------------------------------%
        function n = computeNumSamples(mlmc, tolerance, opt)
            stat = mlmc.getLevelStatistics();
            n = tolerance.^(-2).*sum(sqrt(stat.variance.*stat.cost)) ...
                                   .*sqrt(stat.variance./stat.cost);
            n = ceil(n) - stat.numSamples;
            n = min(n, opt.maxSamples - mlmc.numSamples);
            n = max(n, opt.minSamples - mlmc.numSamples);
            n = min(n, opt.batchSize);
            n = max(n, 0);
        end
        
        %-----------------------------------------------------------------%
        function n = numLevels(mlmc)
            n = numel(mlmc.levels);
        end
        
        %-----------------------------------------------------------------%
        function stat = getLevelStatistics(mlmc)
            lstat = cellfun(@(mcl) mcl.getStatistics, mlmc.levels, 'UniformOutput', false);
            stat = struct();
            for fn = {'estimate', 'variance', 'cost', 'numSamples'}
                stat.(fn{1}) = cellfun(@(lstat) lstat.(fn{1}), lstat);
            end
        end
        
        %-----------------------------------------------------------------%
        % Print
        %-----------------------------------------------------------------%
        function printIterationReport(mlmc, iteration, tolerance, n, rangeStat)
            printIterationReport@MCSimulator(mlmc, iteration, tolerance, n, rangeStat);
            mlmc.printLevelReport(rangeStat);
        end
        
        %-----------------------------------------------------------------%
        function printLevelReport(mlmc, rangeStat)
            header = {'# Samples'     , ...
                      'Estimate'      , ...
                      'Variance'      , ...
                      'Cost'          };
            header = [header; repmat({'It'}, 1, numel(header))];
            header = ['Level'; header(:)]';
            % Field width
            w = max(cellfun(@numel, header), 9);
            % Horizontal line over/under headings and at the end
            hline = @() fprintf([repmat('=', 1, sum(w) + 3*numel(w) + 2), '\n']);
            % Format for each report field
            format = {'u', 'u', 'u', '.2e', '.2e', '.2e', '.2e', '.2e', '.2e'};
            format = cellfun(@(f,n) ['%' num2str(n), f],  ...
                              format, num2cell(w), 'UniformOutput', false);
            % Print header
            fprintf('\n'); hline();
            for i = 1:numel(header)
                fprintf(['| %-' num2str(w(i)), 's '], header{i});
            end
            fprintf(' |\n'); hline();
            for i = 1:mlmc.numLevels()
                % Print field values
                values = struct2cell(mlmc.levels{i}.getStatistics());
                values = values(1:end-1);
                valuesIt = struct2cell(rangeStat{i});
                values   = reshape([values, valuesIt]', [], 1);
                values = [i; values]; %#ok
                for j = 1:numel(header)
                    fprintf(['| ', format{j}, ' '], values{j});
                end
                fprintf(' |\n');
            end
            hline();
        end
        
        %-----------------------------------------------------------------%
        function plotLevels(mlmc, varargin)
            opt = struct('seed', 1);
            [opt, extra] = merge_options(opt, varargin{:});
            bproblem = mlmc.getBaseProblem();
            sample = mlmc.samples.getSample(opt.seed, bproblem);
            for i = 1:mlmc.numLevels
                level   = mlmc.levels{i}.levels{end};
                problem = level.getSampleProblem(opt.seed, 'sample', sample);
                level.setup.plot(problem.SimulatorSetup.model.rock, extra{:});
                plotGrid(level.setup.model.G, 'faceColor', 'none');
                colormap(pink);
            end
        end 
        
    end
    
    methods (Access = protected)
        %-----------------------------------------------------------------%
        function progress = getEnsembleMemberProgress(ensemble, range)
            % Utility function for monitoring the progression of an
            % ensemble member that is being run right now.
            if nargin < 2, range = ensemble.num; end
            progress = zeros(numel(range),1);
            nsteps   = numel(ensemble.setup.schedule.step.val);
            for i = 1:numel(range)
                if exist(fullfile(ensemble.directory(), ...
                        [ensemble.qoi.ResultHandler.dataPrefix, num2str(range(i)), '.mat']), 'file')
                    progress(i) = inf;
                    continue
                end
                dataDir = fullfile(ensemble.directory(), num2str(range(i)));
                if ~exist(dataDir, 'dir')
                    continue;
                end
                files = ls(dataDir);
                progress(i) = numel(folderRegexp(files, 'state\d+\.mat', 'match'))/nsteps;
            end
        end
        
    end
   
end