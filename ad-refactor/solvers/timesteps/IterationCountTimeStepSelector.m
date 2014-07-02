classdef IterationCountTimeStepSelector < SimpleTimeStepSelector
    properties
        targetIterationCount
        maxRelativeAdjustment
        minRelativeAdjustment
    end
    methods
        function selector = IterationCountTimeStepSelector(varargin)
            selector = selector@SimpleTimeStepSelector();
            
            selector.maxRelativeAdjustment = 2;
            selector.minRelativeAdjustment = .5;
            
            selector.targetIterationCount = 10;
            selector = merge_options(selector, varargin{:});
        end
        
        function dt = computeTimestep(selector, dt)
            h = selector.history;
            historyLength = numel(h);
            dt0 = dt;
            dt_new = dt;
            
            if historyLength > 0
                if historyLength > 1 && ...
                   h(end).Iterations > selector.targetIterationCount
                    % We have a bunch of datapoints for iteration counts,
                    % use this for a finite difference estimation
                    dt_1 = h(end).Timestep;
                    dt_0 = h(end-1).Timestep;

                    F_1 = h(end).Iterations;
                    F_0 = h(end-1).Iterations;

                    % Estimate derivative of iterations function based on
                    % timestep length
                    dFdt = (F_1 - F_0)/(dt_1 - dt_0);

                    % Calculate new timestep using this linearization
                    dt_new = dt_1 + (selector.targetIterationCount - F_1)./dFdt;
                end

                if historyLength == 1 || isinf(dt_new) || dt_new < 0;
                    % We either do not have enough data to estimate
                    % convergence rates or we got a bad timestep from
                    % difference approximation
                    prev = h(end).Iterations;
                    dt_prev = h(end).Timestep;
                    if prev > selector.targetIterationCount
                        dt_new = dt_prev/2;
                    elseif prev < selector.targetIterationCount
                        dt_new = 2*dt_prev;
                    end
                end
            else
                % Do nothing, we have no data to work with. Let's hope the
                % first timestep at least converges, so we can do something
                % on the next pass.
            end
            
            change = dt_new/dt;
            change = min(change, selector.maxRelativeAdjustment);
            change = max(change, selector.minRelativeAdjustment);

            dt_new = dt*change;
            
            % Of course, if this timestep is longer than what the
            % non-linear solver wants, who are we to argue...
            dt = min(dt, dt_new);
            dt = max(dt, 0.001*day);

            if selector.verbose && dt ~= dt0
                if historyLength
                    fprintf('Its was: %d ', h(end).Iterations);
                end
                fprintf('Adjusted timestep by a factor %1.2f. dT: %s -> %s\n',...
                    dt/dt0, formatTimeRange(dt0), formatTimeRange(dt));
            end
        end
    end
end