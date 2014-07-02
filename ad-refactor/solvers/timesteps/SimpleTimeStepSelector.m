classdef SimpleTimeStepSelector < handle
    properties
        history
        maxTimestep
        verbose
    end
        
    
    methods
        function selector = SimpleTimeStepSelector(varargin)
            selector.maxTimestep = inf;
            selector.verbose = mrstVerbose();
            
            selector = merge_options(selector, varargin{:});
        end
        
        function reset(selector)
            selector.history = [];
        end
        
        function storeTimestep(selector, report)
            selector.history = vertcat(selector.history, report);
        end
        
        function dt = pickTimestep(selector, dt, model)
            dt = min(selector.maxTimestep, dt);
            dt = selector.computeTimestep(dt);
        end
        
        function dt = computeTimestep(selector, dt) %#ok
            % Trivial case
        end
    end
end