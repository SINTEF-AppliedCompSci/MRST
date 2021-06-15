classdef FlowStateBuilder
    % The FlowStateBuilder class is used to define the AD-state used for
    % computing the fluxes in a simulation. This is done through two
    % mechanisms:
    %  - Limiting of time-steps based on the choice for flow state
    %  - Creating a hybrid state for fluxes where some or all terms are
    %  taken implicitly.
    properties
        verbose = [];
    end
    
    methods
        function fsb = FlowStateBuilder(varargin)
            fsb = merge_options(fsb, varargin{:});
            if isempty(fsb.verbose)
                fsb.verbose = mrstVerbose();
            end
        end

        function dt = getMaximumTimestep(fsb, fd, model, state, state0, dt, forces)
            % Maximum time-step for discretization (required for stability
            % when dealing with explicit methods). The default is implicit
            % in all terms, so we let dt be inf.
            dt = inf;
        end
        
        function flowState = build(builder, fd, model, state, state0, dt)
            % Given two states (representing implicit and explicit
            % choices), hybridize to get a single state where parts are
            % implicit or explicit.
            flowState = state;
        end
        
        function [builder, state] = prepareTimestep(builder, fd, model, state, state0, dt, drivingForces)
            % Called before each nonlinear loop
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
