classdef FactorTimeStepSelector < SimpleTimeStepSelector
    % Time step selector that always tries to increase the time-step
    %
    % SYNOPSIS:
    %   selector = FactorTimeStepSelector();
    %   selector = FactorTimeStepSelector('minRelativeAdjustment', 0.5);
    %   selector = FactorTimeStepSelector('maxRelativeAdjustment', 4);
    %
    % DESCRIPTION:
    %   A simple time step selection class that always tries to multiply
    %   successful time-steps with a factor (default 2) and otherwise cuts
    %   on failure.
    %
    % REQUIRED PARAMETERS:
    %   None
    %
    % OPTIONAL PARAMETERS:
    %   See properties
    %
    % RETURNS:
    %   Selector suitable for passing to the NonLinearSolver class.
    %
    %
    % SEE ALSO:
    %   IterationCountTimeStepSelector, NonLinearSolver, 
    %   StateChangeTimeStepSelector, SimpleTimeStepSelector

    properties

    end
        
    
    methods
        function selector = FactorTimeStepSelector(varargin)
            selector@SimpleTimeStepSelector(varargin{:});
        end

        function dt = cutTimestep(selector, dt_prev, dt, model, solver, state_prev, state_curr, forces)
            dt = dt*selector.minRelativeAdjustment;
        end
        
        function dt = computeTimestep(selector, dt, dt_prev, model, solver, state_prev, state_curr, forces) %#ok
            dt = dt*selector.maxRelativeAdjustment;
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

