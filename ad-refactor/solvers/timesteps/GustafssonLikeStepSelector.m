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
            dt_suggested = dt;
            
            hist = selector.history;
            nHist = numel(hist);
            
            if nHist == 0
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
            
            if hist(end).Converged
                if restart
                    dt_new = (tol/le1)*dt;
                else
                    le0 = hist(end-1).Iterations/solver.maxIterations;
                    dt0 = hist(end-1).Timestep;
                    dt_new = (dt1/dt0)*(tol*le0/le1^2)*dt; 
                end
            else
                if nHist > 1 && ~hist(end-1).Converged
                    % Repeated rejections, bring out the big guns. Probably
                    % not correct for the simple "error" estimation used
                    % here based on iterations
                    le0 = hist(end-1).Iterations/solver.maxIterations;
                    dt0 = hist(end-1).Timestep;
                    p = log(le1/le0)/log(dt1/dt0);
                    dt_new = (tol/le1)^(1/p)*dt;
                else
                    dt_new = (tol/le1)*dt;
                end
            end
            
            dt = dt_new;
            
            if selector.verbose && dt_suggested ~= dt
                fprintf('Prev # its: %d -> ', hist(end).Iterations)
                fprintf('Adjusted timestep by a factor %1.2f. dT: %s -> %s\n',...
                    dt/dt_suggested, formatTimeRange(dt_suggested), formatTimeRange(dt));
            end
        end
    end
end