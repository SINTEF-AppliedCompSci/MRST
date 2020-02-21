classdef ExplicitFlowStateBuilderDG < ExplicitFlowStateBuilder & FlowStateBuilderDG
    methods
        function flowState = build(builder, fd, model, state, state0, dt)
            % Hybridize state
            % The base state is the implicit. Other functions are then
            % assigned.
            flowState = build@FlowStateBuilderDG(builder, fd, model, state, state0, dt);
            flowState = build@ExplicitFlowStateBuilderDG(builder, fd, model, flowState, state0.faceState, dt);
        end
    end
end
