classdef FlowStateBuilder
    % The FlowStateBuilder class is used to define the AD-state used for
    % computing the fluxes in a simulation. This is done through two
    % mechanisms:
    %  - Limiting of time-steps based on the choice for flow state
    %  - Creating a hybrid state for fluxes where some or all terms are
    %  taken implicitly.
    properties
        verbose = mrstVerbose();
    end
    
    methods
        function fsb = FlowStateBuilder(varargin)
            fsb = merge_options(fsb, varargin{:});
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
