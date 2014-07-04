classdef GustafssonLikeStepSelector < SimpleTimeStepSelector
    properties
        targetIterationCount
    end
    methods
        function selector = GustafssonLikeStepSelector(varargin)
            selector = selector@SimpleTimeStepSelector();
            
            selector.targetIterationCount = 5;
                        
            selector = merge_options(selector, varargin{:});
        end
        
        function dt = computeTimestep(selector, dt, model, solver)
            hist = selector.history;
            nHist = numel(hist);
            
            if nHist == 0 || ~hist(end).Converged
                return
            end
            
            
            
            if nHist > 1
                restart = ~hist(end-1).Converged && hist(end).Converged;
            else
                restart = true;
            end
            restart = restart | selector.stepLimitedByHardLimits;
            
            tol = selector.targetIterationCount/solver.maxIterations;

            le1 = hist(end).Iterations/solver.maxIterations;
            dt1 = hist(end).Timestep;
            
            if restart
                dt_new = (tol/le1)*dt;
            else
                le0 = hist(end-1).Iterations/solver.maxIterations;
                dt0 = hist(end-1).Timestep;
                dt_new = (dt1/dt0)*(tol*le0/le1^2)*dt; 
            end
            
            dt = dt_new;
            
        end
    end
end