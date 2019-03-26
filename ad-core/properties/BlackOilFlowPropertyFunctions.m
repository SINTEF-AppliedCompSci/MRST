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
            props.RsMax = RsMax(model, pvt);
            props.RvMax = RvMax(model, pvt);
        end
    end
end
