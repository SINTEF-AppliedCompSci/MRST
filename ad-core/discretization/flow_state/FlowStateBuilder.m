classdef FlowStateBuilder
    properties
        
    end
    
    methods
        function dt = getMaximumTimestep(fsb, fd, model, state, state0, dt, forces)
            dt = inf;
        end
        
        function flowState = build(builder, fd, model, state, state0, dt)
            flowState = state;
        end
    end
end
