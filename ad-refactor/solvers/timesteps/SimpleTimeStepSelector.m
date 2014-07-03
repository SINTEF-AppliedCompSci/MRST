classdef SimpleTimeStepSelector < handle
    properties
        history
        maxTimestep
        minTimestep
        verbose
        maxHistoryLength
        isStartOfCtrlStep
        
        maxRelativeAdjustment
        minRelativeAdjustment
        
        stepLimitedByHardLimits
    end
        
    
    methods
        function selector = SimpleTimeStepSelector(varargin)
            selector.maxTimestep = inf;
            selector.minTimestep = 0;
            selector.maxHistoryLength = 50;
            
            selector.maxRelativeAdjustment = 2;
            selector.minRelativeAdjustment = .5;
            
            selector.verbose = mrstVerbose();
            
            selector = merge_options(selector, varargin{:});
            
            selector.isStartOfCtrlStep = true;
            selector.stepLimitedByHardLimits = true;
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
        
        function dt = pickTimestep(selector, dt, model, solver)
            dt_new = selector.computeTimestep(dt, model, solver);
            
            % Ensure that step does not change too much
            change = dt_new/dt;
            change = min(change, selector.maxRelativeAdjustment);
            change = max(change, selector.minRelativeAdjustment);

            dt = dt*change;
            
            % Honor limits if selection is far off
            dt = min(selector.maxTimestep, dt);
            dt = max(selector.minTimestep, dt);
            
            if dt ~= dt_new
                selector.stepLimitedByHardLimits = true;
            else
                selector.stepLimitedByHardLimits = false;
            end
            selector.isStartOfCtrlStep = false;
        end
        
        function newControlStep(selector, control)
            selector.isStartOfCtrlStep = true;
        end
        
        function dt = computeTimestep(selector, dt, model, solver) %#ok
            % Trivial case
        end
    end
end