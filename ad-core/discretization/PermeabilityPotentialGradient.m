classdef PermeabilityPotentialGradient < GridProperty
    properties
        PermeabilityGradientDiscretization
    end
    
    methods
        function pp = PermeabilityPotentialGradient(backend, kgrad)
            pp@GridProperty(backend);
            pp.PermeabilityGradientDiscretization = kgrad;
        end
        
        function v = evaluateOnGrid(prop, model, state)
            pot = model.getProps(state, 'PhasePotentialDifference');
            nph = numel(pot);
            kgrad = prop.PermeabilityGradientDiscretization;
            v = cell(1, nph);
            for i = 1:nph
                v{i} = kgrad.getPermeabilityGradient(state, pot{i});
            end
        end
    end
end