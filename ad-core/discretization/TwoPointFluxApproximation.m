classdef TwoPointFluxApproximation < PermeabilityGradientDiscretization
    properties
        T
        hasMultiplier
    end
    
    methods
        function tpfa = TwoPointFluxApproximation(model)
            tpfa.T = model.operators.T;
            tpfa.hasMultiplier = isfield(model.fluid, 'pvMult');
        end
        
        function v = getPermeabilityGradient(tpfa, state, value)
            v = tpfa.T.*value;
        end
    end
end