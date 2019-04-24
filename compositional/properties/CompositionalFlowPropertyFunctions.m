classdef CompositionalFlowPropertyFunctions < FlowPropertyFunctions
    properties
        PhaseMixingCoefficients
        ComponentPhaseMassFractions
        ComponentPhaseMoleFractions
        PhaseCompressibilityFactors
        Fugacity
    end
    
    methods
        function props = CompositionalFlowPropertyFunctions(model)
            props@FlowPropertyFunctions(model);
            pvt = props.getRegionPVT(model);
            props.ShrinkageFactors = DensityDerivedShrinkageFactors(model, pvt);
            props.Density = CompositionalDensity(model, pvt);
            
            props.PhaseMixingCoefficients = PhaseMixingCoefficientsLV(model);
            props.Fugacity = FugacityLV(model);
            props.PhaseCompressibilityFactors = PhaseCompressibilityFactorsLV(model);
            props.ComponentPhaseMassFractions = ComponentPhaseMassFractionsLV(model);
            props.ComponentPhaseMoleFractions = ComponentPhaseMoleFractionsLV(model);
            props.Viscosity = CompositionalViscosityLV(model);
        end
    end
end
