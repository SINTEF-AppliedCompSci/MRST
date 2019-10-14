classdef FlowStateBuilderDG < FlowStateBuilder
    
    methods
        function flowState = build(builder, fd, model, state, state0, dt)
            
            flowState = state.faceStateDG;
            propfn = model.getStateFunctionGroupings();
            for i = 1:numel(propfn)
                p = propfn{i};
                struct_name = p.getStateFunctionContainerName();
                if isfield(state, struct_name)
                    flowState = rmfield(flowState, struct_name);
                end
            end
            flowState = model.initStateFunctionContainers(flowState);
%             state = model.validateState(state);
            flowState.type = 'face';
            if ~isfield(flowState, 'cells')
                flowState.cells = (1:model.G.cells.num)';
                flowState.faces = (1:model.G.faces.num)';
            end
        end
    end
end