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
                % Level QoIs will be stored in a subfolder `level-<levelNo>`
                directory = fullfile(mlmc.getDataPath(), ...
                            ['level-', num2str(levels{lix(end)}.levelNo)]);
                mlmc.levels{i} = MCLevelSimulator(setup, samples, qoi, levels(lix), ...
                                                  varargin{:}           , ...
                                                  'directory', directory);
            end
        end
        
        %-----------------------------------------------------------------%
        function range = runBatch(mlmc, varargin)
            opt = struct('batchSize', []);
            [opt, extra] = merge_options(opt, varargin{:});
            range = nan(sum(opt.batchSize), 1);
            if mlmc.verbose
                if mlmc.numSamples == 0
                    s = 'warmup';
                else
                    s = 'iteration';
                end
                bs = num2str(opt.batchSize');
                fprintf(['Running multilevel Monte Carlo %s ', ...
                         'with (%s) samples\n\n'], s, bs     );
            end
            for i = 1:mlmc.numLevels
                if opt.batchSize(i) > 0
                    maxId = mlmc.getLargestSeed();
                    lrange = (1:opt.batchSize(i)) + maxId;
                    mlmc.levels{i}.runBatch('range', lrange, extra{:});
                    range((1:opt.batchSize(i)) + sum(opt.batchSize(1:i-1))) = lrange;
                    while ~all(mlmc.levels{i}.getSimulationStatus(lrange) > 0)
                        pause(0.05);
                    end
                    rangeStat = mlmc.levels{i}.updateStatistics(lrange);
                else
                    rangeStat = [];
                end
                mlmc.printLevelReport(i, rangeStat);
            end
        end
        
        %-----------------------------------------------------------------%
        function rangeStat = updateStatistics(mlmc, range)
            rangeStat = cell(mlmc.numLevels,1);
            if ~isempty(range)
                for i = 1:mlmc.numLevels()
                    ids = mlmc.levels{i}.qoi.ResultHandler.getValidIds();
                    lrange = range(ismember(range, ids));
                    if ~isempty(lrange)
                        rangeStat{i} = mlmc.levels{i}.updateStatistics(lrange);
                    end
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
        function reset(mlmc, varargin)
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
                    return
                end
            end
            for i = 1:mlmc.numLevels
                mlmc.levels{i}.reset('prompt', false, extra{:});
            end
            reset@MCSimulator(mlmc, 'prompt', false);
        end
        
        %-----------------------------------------------------------------%
        function resetStatistics(mlmc)
            resetStatistics@MCSimulator(mlmc);
            cellfun(@(level) level.resetStatistics(), mlmc.levels);
        end
        
        %-----------------------------------------------------------------%
        function [n, n0] = computeNumSamples(mlmc, tolerance, opt)
            stat = mlmc.getLevelStatistics();
            n0 = tolerance.^(-2).*sum(sqrt(stat.variance.*stat.cost)) ...
                                    .*sqrt(stat.variance./stat.cost);
            n0 = max(ceil(n0) - stat.numSamples, 0);
            n  = min(n0, opt.maxSamples - mlmc.numSamples);
            n  = max(n , opt.minSamples - mlmc.numSamples);
            n  = min(n , opt.batchSize);
            n  = max(n , 0);
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
        function printIterationReport(mlmc, iteration, tolerance, n)
            printIterationReport@MCSimulator(mlmc, iteration, tolerance, 0, true);
            
        end
        
        %-----------------------------------------------------------------%
        function printLevelReport(mlmc, levelNo, rangeStat)
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
            if levelNo == 1
                % Print header
                hline();
                for i = 1:numel(header)
                    fprintf(['| %-' num2str(w(i)), 's '], header{i});
                end
                fprintf(' |\n'); hline();
            end
            % Print field values
            values = struct2cell(mlmc.levels{levelNo}.getStatistics());
            values = values(1:end-1);
            if ~isempty(rangeStat)
                valuesIt = struct2cell(rangeStat);
            else
                valuesIt = num2cell(nan(size(values)));
            end
            values   = reshape([values, valuesIt]', [], 1);
            values = [levelNo; values];
            for j = 1:numel(header)
                fprintf(['| ', format{j}, ' '], values{j});
            end
            fprintf(' |\n');
            if levelNo == mlmc.numLevels, hline(); end
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
