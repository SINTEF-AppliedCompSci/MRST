classdef BlackOilFlowPropertyFunctions < FlowPropertyFunctions
    properties
        RsMax
        RvMax
    end
    
    methods
        function props = BlackOilFlowPropertyFunctions(model)
            props@FlowPropertyFunctions(model);
            % Get PVT region from density
            pvt = props.Density.regions;
            props.RsMax = RsMax(model.AutoDiffBackend, pvt);
            props.RvMax = RvMax(model.AutoDiffBackend, pvt);
        end
        
        function state = evaluateProperty(props, model, state, name)
            switch name
                case 'RvMax'
                    state = props.evaluateDependencies(model, state, {'CapillaryPressure'});
            end
            state = evaluateProperty@FlowPropertyFunctions(props, model, state, name);
        end
    end
end
