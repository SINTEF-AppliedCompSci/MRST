classdef SimpleTimeStepSelector < handle
    % Time step selector base class
    %
    % SYNOPSIS:
    %   selector = SimpleTimeStepSelector();
    %   selector = SimpleTimeStepSelector('maxTimestep', 5*day);
    %
    % DESCRIPTION:
    %   The timestep selector base class is called by the NonLinearSolver to
    %   determine timesteps, based on hard limits such as the min/max timesteps
    %   as well as possibly more advanced features via the computeTimestep
    %   method that can account for iteration count, residual reduction etc.
    %
    % REQUIRED PARAMETERS:
    %   None
    %
    % OPTIONAL PARAMETERS:
    %   See properties
    %
    % RETURNS:
    %   selector suitable for passing to the NonLinearSolver class.
    %
    %
    % SEE ALSO:
    %   IterationCountTimeStepSelector, NonLinearSolver, 
    %   StateChangeTimeStepSelector

    properties
        history % Stored history used to pick next timestep
        maxTimestep % Hard upper limit on timestep in seconds
        minTimestep % Hard lower limit on timestep in seconds
        verbose % Extra output
        maxHistoryLength % The maximum number of history steps stored
        isStartOfCtrlStep % Flag indicating the beginning of a control step
        controlsChanged % Flag indicating that the controls have changed
        
        firstRampupStep % The first ministep attempted after controls have changed
        % Could be set to a low value to get the timestep controller started 
        % with some estimate of problem stiffness.
        % 
        firstRampupStepRelative % Relative version of firstRampupStep
        % (i.e. if timestep is 5*days, the ramp up step will then be
        % 5*days*firstRampupStepRelative.

        maxRelativeAdjustment % Ensure dt_next < dt_suggested*maxRelativeAdjustment
        minRelativeAdjustment % Ensure dt_next > dt_suggested*minRelativeAdjustment
        
        stepLimitedByHardLimits % Flag indicating that hard limits and not any step algorithm was the cause of the previous timestep taken
        previousControl % Previous control seen by the selector used to determine when controls change.
        resetOnControlsChanged % Reset when controls change
        isFirstStep % Flag indicating that we are at the start of the simulation
    end
        
    
    methods
        function selector = SimpleTimeStepSelector(varargin)
            selector.maxTimestep = inf;
            selector.minTimestep = 0;
            selector.maxHistoryLength = 50;
            
            selector.maxRelativeAdjustment = 2;
            selector.minRelativeAdjustment = .5;
            
            selector.firstRampupStep = inf;
            selector.firstRampupStepRelative = 1;
            
            selector.resetOnControlsChanged = false;

            selector.verbose = mrstVerbose();
            
            selector = merge_options(selector, varargin{:});
            
            assert(selector.firstRampupStepRelative <= 1 && ...
                   selector.firstRampupStepRelative >  0);
            
            selector.reset();
        end
        
        function reset(selector)
            selector.history = [];
            selector.previousControl = [];
            selector.isStartOfCtrlStep = true;
            selector.isFirstStep       = true;
            selector.controlsChanged   = true;
            selector.stepLimitedByHardLimits = true;
        end
        
        function storeTimestep(selector, report)
            if isempty(selector.history)
                selector.history = report;
            else
                selector.history = vertcat(selector.history, report);
            end
            
            n = selector.maxHistoryLength;
            if numel(selector.history) > n
                selector.history = selector.history((end-n):end);
            end
        end
        
        function dt = pickTimestep(selector, dt_prev, dt, model, solver, state_prev, state_curr, forces)
            if selector.controlsChanged && ...
              (selector.resetOnControlsChanged || selector.isFirstStep)
                % First relative check
                dt = min(dt, selector.firstRampupStepRelative*dt);
                % Then absolute check
                dt = min(dt, selector.firstRampupStep);
                selector.stepLimitedByHardLimits = true;
            end
            dt0 = dt;
            dt_new = selector.computeTimestep(dt, dt_prev, model, solver, state_prev, state_curr, forces);

            % Ensure that step does not change too much
            change = dt_new/dt;
            if ~selector.isStartOfCtrlStep
                % Apply maximum relative adjustment only if we are not at
                % start of control step
                change = min(change, selector.maxRelativeAdjustment);
                change = max(change, selector.minRelativeAdjustment);
            end

            dt = dt*change;
            
            % Honor limits if selection is far off
            dt = min(selector.maxTimestep, dt);
            dt = max(selector.minTimestep, dt);
            
            if selector.verbose && dt0 ~= dt
                if ~isempty(selector.history)
                    fprintf('Prev # its: %d -> ', selector.history(end).Iterations)
                end
                fprintf('Suggested timestep modification of a factor %1.2f. dT: %s -> %s\n',...
                    dt/dt0, formatTimeRange(dt0, 2), formatTimeRange(dt, 2));
            end

            if dt ~= dt_new
                selector.stepLimitedByHardLimits = true;
            else
                selector.stepLimitedByHardLimits = false;
            end
            % If we originally were at the start of a control step, we are
            % no longer at the beginning. Reset those indicators.
            selector.isStartOfCtrlStep = false;
            selector.isFirstStep = false;
            selector.controlsChanged = false;
        end
        
        function newControlStep(selector, control)
            % Determine if controls have changed
            selector.isStartOfCtrlStep = true;
            
            prev = selector.previousControl;
            
            if isempty(prev) || prev.controlId ~= control.controlId
                % We are at the beginning of a new control step and should
                % reset the selector accordingly
                if selector.resetOnControlsChanged
                    selector.reset();
                end
                selector.controlsChanged = true;
                selector.previousControl = control;
            else
                selector.controlsChanged = false;
            end
        end

        function dt = cutTimestep(selector, dt_prev, dt, model, solver, state_prev, state_curr, forces)
            dt = dt*selector.minRelativeAdjustment;
        end
        
        function dt = computeTimestep(selector, dt, dt_prev, model, solver, state_prev, state_curr, forces) %#ok
            % Compute timestep dynamically - does nothing for base class    
        end
    end
end

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

