classdef CompositionalFlowPropertyFunctions < FlowPropertyFunctions
    properties
        ComponentPhaseMassFractions
        ComponentPhaseMoleFractions
    end
    
    methods
        function props = CompositionalFlowPropertyFunctions(model)
            props@FlowPropertyFunctions(model);
            pvt = props.getRegionPVT(model);
            props = props.setStateFunction('ShrinkageFactors', DensityDerivedShrinkageFactors(model, pvt));
            props = props.setStateFunction('Density', CompositionalDensity(model, pvt));
            
            props = props.setStateFunction('ComponentPhaseMassFractions', ComponentPhaseMassFractionsLV(model));
            props = props.setStateFunction('ComponentPhaseMoleFractions', ComponentPhaseMoleFractionsLV(model));

            props = props.setStateFunction('PhaseMixingCoefficients', PhaseMixingCoefficientsLV(model));
            props = props.setStateFunction('Fugacity', FugacityLV(model));
            props = props.setStateFunction('PhaseCompressibilityFactors', PhaseCompressibilityFactorsLV(model));
            props = props.setStateFunction('Viscosity', CompositionalViscosityLV(model));
        end
        
        function props = setCompactEvaluation(props, val)
            for i = 1:numel(props.propertyNames)
                name = props.propertyNames{i};
                fn = props.getStateFunction(name);
                if isprop(fn, 'useCompactEvaluation')
                    fn.useCompactEvaluation = val;
                    props = props.setStateFunction(name, fn);
                end
            end
        end
    end
end
