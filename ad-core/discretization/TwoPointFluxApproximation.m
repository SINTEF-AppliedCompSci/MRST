classdef TwoPointFluxApproximation < PermeabilityGradientDiscretization
    properties
        T
    end
    
    methods
        function tpfa = TwoPointFluxApproximation(model)
            tpfa.T = model.operators.T;
        end
        
        function v = getPermeabilityGradient(tpfa, state, value)
            v = tpfa.T.*value;
        end
    end
end