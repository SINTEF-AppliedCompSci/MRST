classdef ExplicitFlowStateBuilderDG < ExplicitFlowStateBuilder & FlowStateBuilderDG
    methods
        function flowState = build(builder, fd, model, state, state0, dt, type)
            % Hybridize state
            % The base state is the implicit. Other functions are then
            % assigned.
            flowState = build@FlowStateBuilderDG(builder, fd, model, state, state0, dt, type);
            flowState = build@ExplicitFlowStateBuilder(builder, fd, model, flowState, state0.faceStateDG, dt);
        end
    end
end
