classdef PermeabilityPotentialGradient < GridProperty
    properties
        PermeabilityGradientDiscretization
    end
    
    methods
        function pp = PermeabilityPotentialGradient(backend, kgrad)
            pp@GridProperty(backend);
            pp.PermeabilityGradientDiscretization = kgrad;
            pp = pp.dependsOn('PhasePotentialDifference');
            if isa(kgrad, 'TwoPointFluxApproximation')
                pp = pp.dependsOn('Transmissibility');
            end
        end
        
        function v = evaluateOnDomain(prop, model, state)
            pot = prop.getEvaluatedDependencies(state, 'PhasePotentialDifference');
            nph = numel(pot);
            kgrad = prop.PermeabilityGradientDiscretization;
            v = cell(1, nph);
            for i = 1:nph
                v{i} = kgrad.getPermeabilityGradient(model, state, pot{i});
            end
        end
    end
end