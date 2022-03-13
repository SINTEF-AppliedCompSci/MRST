classdef MCSimulator < MRSTEnsemble
    
    properties(Access = protected)
        numSamples = 0
        estimate   = 0
        variance   = 0
        cost       = 0
        rmse       = inf
        history    = {}
        included
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function mc = MCSimulator(mrstExample, samples, qoi, varargin)
            assert(numel(qoi.names) == 1, 'MCSimulator only supports one qoi');
            mc = mc@MRSTEnsemble(mrstExample, samples, qoi, varargin{:});
            mc.figures.progressMC = [];
        end
        
        %-----------------------------------------------------------------%
        function runSimulation(mc, varargin)
            opt = struct('range'       , inf , ...
                         'batchSize'   , inf , ...
                         'maxSamples'  , 10  , ...
                         'minSamples'  , 2   , ...
                         'tolerance'   , -inf, ...
                         'relTolerance', -inf, ...
                         'prompt'      , true, ...
                         'plotProgress', true);
            opt = merge_options(opt, varargin{:});
            assert(xor(isfinite(opt.tolerance), isfinite(opt.relTolerance)), ...
                   'Choose either relative or absolute tolerance');
            % Check for computed samples
            range = mc.getComputedSampleRange();
            n     = numel(range);
            if n > 0 && opt.prompt
                % We found computed samples - ask user if they should be
                % included in the estimate
                prompt = sprintf(['Found %d QoIs for this problem. '   , ...
                                  'Would you like to include these in ', ...
                                  'your estimate (y)? If not, I will ' , ...
                                  'delete them y/n [y]: '], n          );
                str = input(prompt,'s');
                if strcmpi(str, 'y') || isempty(str)
                    % Include samples in estimate
                    fprintf('Ok, I will include the QoIs in the estimate.\n');
                    if mc.numSamples ~= n
                        mc.resetStatistics();
                        mc.updateStatistics(range);
                    end
                else
                    % Reset simulator. This includes wiping out any samples
                    % already computed
                    mc.reset('prompt', false);
                end
            end
            % Set tolerance
            tolerance = opt.tolerance;
            if isinf(opt.tolerance)
                % Set tolerance based on relative reduction target
                tolerance = opt.relTolerance*mc.estimate;
            end
            % Compute number of samples needed
            n = mc.computeNumSamples(tolerance, opt);
            iteration = 1;
            while (any(n > 0) && mc.numSamples < opt.maxSamples) ... 
                              || mc.numSamples < opt.minSamples
                % Run batch of n new samples
                range = mc.runBatch('batchSize', n, 'plotProgress', opt.plotProgress);
                % Update statistics
                mc.updateStatistics(range);
                if isinf(opt.tolerance)
                    % Adjust tolerance if we aim at relative reduction
                    tolerance = opt.relTolerance*mc.estimate;
                end
                % Recompute number of samples needed
                [n, n0] = mc.computeNumSamples(tolerance, opt);
                if mc.verbose
                    mc.printIterationReport(iteration, tolerance, n)
                    if any(n < n0)
                        if numel(n) == 1
                            s = 'maximum is';
                        else
                            s = 'maxima are';
                        end
                        warning(['Requested (%s) samples, but '  , ...
                                 'prescribed %s (%s). Reducing ' , ...
                                 'number of samples'            ], ...
                                 num2str(n0'), s, num2str(opt.batchSize'))
                    end
                end
                if opt.plotProgress
                    % Plot progress of estimate with rmse bounds
                    out = mc.getHistory();
                    mc.figures.progressMC = plotMonteCarloProgress(out, ...
                                                    mc.figures.progressMC);
                    drawnow(); pause(0.1);
                end
                % Update iteration count
                iteration = iteration + 1;
            end
            
        end
        
        %-----------------------------------------------------------------%
        function range = runBatch(mc, varargin)
            range = mc.simulateEnsembleMembers(varargin{:});
        end
        
        %-----------------------------------------------------------------%
        function maxId = getLargestSeed(mc)
            maxId = max(mc.qoi.ResultHandler.getValidIds());
            if isempty(maxId), maxId = 0; end
        end
        
        %-----------------------------------------------------------------%
        function range = getComputedSampleRange(mc)
            range = mc.qoi.ResultHandler.getValidIds();
        end
        
        %-----------------------------------------------------------------%
        function reset(mc, varargin)
            reset@MRSTEnsemble(mc, varargin{:});
            mc.resetStatistics();
        end
        
        %-----------------------------------------------------------------%
        function resetStatistics(mc)
            mc.estimate   = 0;
            mc.variance   = 0;
            mc.cost       = 0;
            mc.numSamples = 0;
            mc.rmse       = inf;
            mc.history    = {};
            mc.included   = [];
        end
        
        %-----------------------------------------------------------------%
        function rangeStat = updateStatistics(mc, range)
            if nargin < 2 || all(isinf(range))
                range = mc.qoi.ResultHandler.getValidIds();
                mc.resetStatistics();
            end
            range = range(~ismember(range, mc.included));
            range = range(mc.getSimulationStatus(range)>0);
            ids = inf;
            while ~all(ismember(range, ids))
                ids = mc.qoi.ResultHandler.getValidIds();
                pause(0.05);
            end
            if isempty(range), rangeStat = []; return; end
            % Get current statistics
            m0 = mc.estimate;
            v0 = mc.variance;
            c0 = mc.cost;
            n0 = mc.numSamples;
            % Get batch statistics
            [m, v] = mc.qoi.computeQoIMean(range);
            c = m.cost; m = m.(mc.qoi.names{1}); v = v.(mc.qoi.names{1});
            n = numel(range);
            % Return statistics for this range
            rangeStat = struct('numSamples', n, ...
                               'estimate'  , m, ...
                               'variance'  , v, ...
                               'cost'      , c);
            % Update statistics
            mc.estimate   = computeMean(m0, m, n0, n);
            mc.variance   = computeVariance(v0, v, m0, m, n0, n);
            mc.cost       = computeMean(c0, c, n0, n);
            mc.numSamples = mc.numSamples + numel(range);
            if mc.numSamples > 1
                mc.rmse = sqrt(mc.variance/mc.numSamples);
            end
            mc.included = [mc.included, reshape(range, 1, [])];
            % Update history
            mc.history{end+1} = struct('numSamples', mc.numSamples, ...
                                       'estimate'  , mc.estimate  , ...
                                       'variance'  , mc.variance  , ...
                                       'cost'      , mc.cost      , ...
                                       'rmse'      , mc.rmse      );
        end
        
        %-----------------------------------------------------------------%
        function [n, n0] = computeNumSamples(mc, tolerance, opt)
            n0 = ceil(mc.variance/tolerance^2) - mc.numSamples;
            n = min(n0, opt.maxSamples - mc.numSamples);
            n = max(n , opt.minSamples - mc.numSamples);
            n = min(n , opt.batchSize);
            n = max(n , 0);
        end
        
        %-----------------------------------------------------------------%
        % Getters                                                         %
        %-----------------------------------------------------------------%
        function stat = getStatistics(mc)
            stat = struct();
            stat.numSamples = mc.numSamples;
            stat.estimate   = mc.estimate;
            stat.variance   = mc.variance;
            stat.cost       = mc.cost;
            stat.rmse       = mc.rmse;
        end
        
        %-----------------------------------------------------------------%
        function history = getHistory(mc, batchSize)
            if nargin > 1
                mc.assembleHistory(batchSize);
            end
            get = @(fn) reshape(cellfun(@(h) h.(fn), mc.history), [], 1);
            history = struct();
            history.numSamples = get('numSamples');
            history.estimate   = get('estimate');
            history.variance   = get('variance');
            history.cost       = get('cost');
            history.rmse       = get('rmse');
        end
        
         %-----------------------------------------------------------------% 
        function assembleHistory(mc, batchSize)
            mc.resetStatistics();
            stat = {nan};
            i = 0;
            while ~all(cellfun(@isempty, stat))
                range = (1:batchSize) + i*batchSize;
                stat = mc.updateStatistics(range);
                if ~iscell(stat), stat = {stat}; end
                i = i + 1;
            end
        end
        
        %-----------------------------------------------------------------%
        % Print
        %-----------------------------------------------------------------%
        function printIterationReport(mc, iteration, tolerance, n, printHeader)
            if nargin < 5, printHeader = false; end
            header = {'Iteration' , ...
                      '# Samples', ...
                      'Estimate'  , ...
                      'Variance'  , ...
                      'Cost'      , ...
                      'RMSE'      , ...
                      'Tolerance' };
            % Field width
            w = max(cellfun(@numel, header), 8);
            % Horizontal line over/under headings and at the end
            hline = @() fprintf([repmat('=', 1, sum(w) + 3*numel(w) + 2), '\n']);
            % Format for each report field
            format = {'u', 'u', '.2e', '.2e', '.2f', '.2e', '.2e'};
            format = cellfun(@(f,n) ['%' num2str(n), f],  ...
                              format, num2cell(w), 'UniformOutput', false);
            if iteration == 1 || printHeader
                % Print header
                fprintf('\n'); hline();
                for i = 1:numel(header)
                    fprintf(['| %-' num2str(w(i)), 's '], header{i});
                end
                fprintf(' |\n'); hline();
            end
            % Print field values
            values = struct2cell(mc.getStatistics);
            values = [iteration; values; tolerance];
            for i = 1:numel(header)
                fprintf(['| ', format{i}, ' '], values{i});
            end
            fprintf(' |\n');
            % Print 
            if ~any(n > 0), hline(); fprintf('\n'); return; end
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
