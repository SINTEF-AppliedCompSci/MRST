classdef SimpleTimeStepSelector < handle
    properties
        history
        maxTimestep
        minTimestep
        verbose
        maxHistoryLength
        isStartOfCtrlStep
    end
        
    
    methods
        function selector = SimpleTimeStepSelector(varargin)
            selector.maxTimestep = inf;
            selector.minTimestep = 0;
            selector.maxHistoryLength = 50;
            
            selector.verbose = mrstVerbose();
            
            selector = merge_options(selector, varargin{:});
            
            selector.controlsChanged = true;
        end
        
        function reset(selector)
            selector.history = [];
        end
        
        function storeTimestep(selector, report)
            selector.history = vertcat(selector.history, report);
            
            n = selector.maxHistoryLength;
            if numel(selector.history) > n
                selector.history = selector.history((end-n):end);
            end
        end
        
        function dt = pickTimestep(selector, dt, model)
            dt = selector.computeTimestep(dt);
            
            % Honor limits if selection is far off
            dt = min(selector.maxTimestep, dt);
            dt = max(selector.minTimestep, dt);
            
            selector.isStartOfCtrlStep = false;
        end
        
        function newControlStep(selector, control)
            selector.isStartOfCtrlStep = true;
        end
        
        function dt = computeTimestep(selector, dt) %#ok
            % Trivial case
        end
    end
end