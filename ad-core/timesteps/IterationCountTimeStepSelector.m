classdef IterationCountTimeStepSelector < SimpleTimeStepSelector
    % Adjust timesteps based with target iteration count, based on history
    %
    % SYNOPSIS:
    %   selector = IterationCountTimeStepSelector();
    %   selector = IterationCountTimeStepSelector('targetIterationCount', 5);
    %
    % DESCRIPTION:
    %   Routine used for dynamic timestepping 
    %
    % REQUIRED PARAMETERS:
    %   None.
    %
    % OPTIONAL PARAMETERS:
    %   targetIterationCount  - Desired number of iterations.
    %
    %   iterationOffset       - Uses [actual + iterationOffset] to calculate
    %                           the parameter. Larger values makes the step
    %                           selector less aggressive for iteration targets
    %                           near zero.
    %
    %   (Other options)       - Inherited from SimpleTimeStepSelector.
    %
    % RETURNS:
    %   Time step selector.
    %
    %
    % SEE ALSO:
    %   SimpleTimeStepSelector, NonLinearSolver

    properties
        targetIterationCount % Desired number of nonlinear iterations per timestep.
        iterationOffset % Offset to make iteration a bit smoother as a response function.
    end
    methods
        function selector = IterationCountTimeStepSelector(varargin)
            selector = selector@SimpleTimeStepSelector();
            
            selector.targetIterationCount = 5;
            selector.iterationOffset = 5;
            
            selector = merge_options(selector, varargin{:});
        end
        
        function dt = computeTimestep(selector, dt, dt_prev, model, solver, state_prev, state_curr, forces)
            % Dynamically compute timestep
            hist = selector.history;
            nHist = numel(hist);
            
            if nHist == 0 || ~hist(end).Converged ||...
               (selector.controlsChanged && selector.resetOnControlsChanged)
                return
            end
            
            if model.stepFunctionIsLinear
                if isa(model, 'SequentialPressureTransportModel')
                    % Use transport solver as the iteration counter as
                    % the pressure equation should be less nonlinear.
                    getIts = @(x) x.NonlinearReport{end}.TransportSolver.Iterations;
                else
                    error(['Step function is linear, but I do not know',...
                        ' how to calculate the iterations for models of type ', ...
                        class(model)]);
                end 
            else
                getIts = @(x) x.Iterations;
            end
            
            if nHist > 1
                % We consider this a restart if we were limited by bounds
                % outside of the function or if the previous step was the
                % first that converged
                restart = selector.stepLimitedByHardLimits |...
                          ~hist(end-1).Converged;
            else
                restart = true;
            end
            maxits = solver.maxIterations + selector.iterationOffset;
            offset = selector.iterationOffset;
            
            tol = (selector.targetIterationCount + offset)/maxits;

            le1 = (getIts(hist(end)) + offset)/maxits;
            dt1 = hist(end).Timestep;
            
            if restart
                dt_new = (tol/le1)*dt1;
            else
                le0 = (getIts(hist(end-1)) + offset)/maxits;
                dt0 = hist(end-1).Timestep;
                dt_new = (dt1/dt0)*(tol*le0/le1^2)*dt1; 
            end
            
            dt = dt_new;
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

