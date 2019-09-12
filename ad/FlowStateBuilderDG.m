classdef FlowStateBuilderDG < FlowStateBuilder
    
    methods
        function flowState = build(builder, fd, model, state, state0, dt)
            propfn = model.getStateFunctionGroupings();
            for i = 1:numel(propfn)
                p = propfn{i};
                struct_name = p.getStateFunctionContainerName();
                if isfield(state, struct_name)
                    state = rmfield(state, struct_name);
                end
            end
%             state = model.validateState(state);
            flowState = state;
            flowState.type = 'face';
        end
    end
end