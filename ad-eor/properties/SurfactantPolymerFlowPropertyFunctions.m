classdef SurfactantPolymerFlowPropertyFunctions < SurfactantFlowPropertyFunctions

    properties
        PolymerAdsorption;
        PolymerViscMultiplier;
        EffectiveMixturePolymerViscMultiplier;
    end
    
    methods
        function props = SurfactantPolymerFlowPropertyFunctions(model)
            props = props@SurfactantFlowPropertyFunctions(model);
            sat = props.getRegionSaturation(model);
            props.PolymerAdsorption     = PolymerAdsorption(model, sat);
            props.PolymerViscMultiplier = PolymerViscMultiplier(model, sat);
            props.EffectiveMixturePolymerViscMultiplier = EffectiveMixturePolymerViscMultiplier(model, sat);
            props.Viscosity             = SurfactantPolymerViscosity(model, sat);
        end
    end
end