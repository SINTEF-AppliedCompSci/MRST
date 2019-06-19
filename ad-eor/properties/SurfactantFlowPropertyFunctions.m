classdef SurfactantFlowPropertyFunctions < FlowPropertyFunctions
    properties
        CapillaryNumber
    end
    
    methods
        function props = SurfactantFlowPropertyFunctions(model)
            props = props@FlowPropertyFunctions(model);
            props.RelativePermeability = SurfactantRelativePermeability(model);
            props.CapillaryPressure = SurfactantCapillaryPressure(model);
            props.CapillaryNumber = CapillaryNumber(model);
        end
    end
end