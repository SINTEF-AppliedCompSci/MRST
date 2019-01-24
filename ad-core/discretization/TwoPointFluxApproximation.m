classdef TwoPointFluxApproximation < PermeabilityGradientDiscretization
    properties
    end
    
    methods
        function tpfa = TwoPointFluxApproximation(model)
            
        end

        function v = getPermeabilityGradient(tpfa, model, state, value)
            T = model.getProp(state, 'transmissibility');
            v = T.*value;
        end
    end
end