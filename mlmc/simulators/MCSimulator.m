classdef MCSimulator < MRSTEnsemble
    properties(Access = protected)
        estimate   = 0
        variance   = 0
        cost       = 0
        numSamples = 0
        rmse       = inf
        history    = {}
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function runMonteCarloSimulation(mc, varargin)
            opt = struct('range'       , inf , ...
                         'batchSize'   , inf , ...
                         'maxSamples'  , 10  , ...
                         'minSamples'  , 2   , ...
                         'tolerance'   , -inf, ...
                         'relTolerance', 1e-2, ...
                         'prompt'      , true, ...
                         'plotProgress', true);
            opt = merge_options(opt, varargin{:});
            
            n = mc.qoi.ResultHandler.numelData();
            
            if n > 0 && opt.prompt
                prompt = sprintf(['Found %d QoIs for this problem. '              , ...
                                  'Would you like to include these in '           , ...
                                  'your estimate (y)? If not, I will delete them ', ...
                                  'y/n [y]: '], n                                 );
                str = input(prompt,'s');
                if strcmpi(str, 'y') || isempty(str)
                    fprintf('Ok, I will include the QoIs in the estimate.\n');
                else
                    mc.reset('prompt', false);
                end
            end
            
            if mc.numSamples ~= mc.qoi.ResultHandler.numelData()
                mc.resetStatistics();
                range = mc.qoi.ResultHandler.getValidIds();
                if ~isempty(range)
                    mc.updateStatistics(range);
                end
            end
            
            tolerance = opt.tolerance;
            if isinf(opt.tolerance)
                tolerance = opt.relTolerance*mc.estimate;
            end
            n = mc.computeNumSamples(tolerance, opt);
            h_mcprogress = [];
            while (n > 0                           && ...
                   mc.numSamples < opt.maxSamples) || ...
                   mc.numSamples < opt.minSamples
                maxId = max(mc.qoi.ResultHandler.getValidIds());
                if isempty(maxId), maxId = 0; end
                range = (1:n) + maxId;
                [h_progress, h_qoi] = mc.simulateEnsembleMembers('range', range, 'plotProgress', opt.plotProgress);
                
                mc.updateStatistics(range);
                if isinf(opt.tolerance)
                    tolerance = opt.relTolerance*mc.estimate;
                end
                n = mc.computeNumSamples(tolerance, opt);
                if opt.plotProgress
                    out = mc.getHistory();
                    delete(h_progress); delete(h_qoi);
                    h_mcprogress = plotMonteCarloProgress(out, h_mcprogress);
                    drawnow(); pause(0.1);
                end
            end
            
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
        end
        
        %-----------------------------------------------------------------%
        function updateStatistics(mc, range)
            if nargin < 2 || all(isinf(range))
                range = mc.qoi.ResultHandler.getValidIds();
                mc.resetStatistics();
            end
            % Get current statistics
            m0 = mc.estimate;
            v0 = mc.variance;
            c0 = mc.cost;
            n0 = mc.numSamples;
            % Get batch statistics
            [m, v] = mc.qoi.computeMean(range); m = m{1};
            c = 1;
            n = numel(range);
            % Update statistics
            mc.estimate   = computeMean(m0, m, n0, n);
            mc.variance   = computeVariance(v0, v, m0, m, n0, n);
            mc.cost       = computeMean(c0, c, n0, n);
            mc.numSamples = mc.numSamples + numel(range);
            if mc.numSamples > 1
                mc.rmse = sqrt(mc.variance/mc.numSamples);
            end
            % Update history
            mc.history{end+1} = struct('estimate'  , mc.estimate  , ...
                                       'variance'  , mc.variance  , ...
                                       'cost'      , mc.cost      , ...
                                       'numSamples', mc.numSamples, ...
                                       'rmse'      , mc.rmse      );
        end
        
        %-----------------------------------------------------------------%
        function n = computeNumSamples(mc, tolerance, opt)
            n = ceil(mc.variance/tolerance^2) - mc.numSamples;
            n = min(n, opt.batchSize);
            n = min(n, opt.maxSamples - mc.numSamples);
            n = max(n, opt.minSamples - mc.numSamples);
            n = max(n, 0);
        end
        
        %-----------------------------------------------------------------%
        % Getters                                                         %
        %-----------------------------------------------------------------%
        function out = getStatistics(mc)
            out = struct();
            out.estimate   = mc.estimate;
            out.variance   = mc.variance;
            out.cost       = mc.cost;
            out.numSamples = mc.numSamples;
            out.rmse       = mc.rmse;
        end
        
        %-----------------------------------------------------------------%
        function out = getHistory(mc)
            get = @(fn) reshape(cellfun(@(h) h.(fn), mc.history), [], 1);
            out = struct();
            out.estimate   = get('estimate');
            out.variance   = get('variance');
            out.cost       = get('cost');
            out.numSamples = get('numSamples');
            out.rmse       = get('rmse');
        end

    end
    
end