classdef BlackOilFlowPropertyFunctions < FlowPropertyFunctions
    properties
        RsMax
        RvMax
    end
    
    methods
        function props = BlackOilFlowPropertyFunctions(model)
            props@FlowPropertyFunctions(model);
        end
        
        function state = evaluateProperty(props, model, state, name)
            switch name
                case 'rvMax'
                    state = props.evaluateDependencies(model, state, {'CapillaryPressure'});
            end
            state = evaluateProperty@FlowPropertyFunctions(props, model, state, name);
        end

    end
end
