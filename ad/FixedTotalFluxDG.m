classdef FixedTotalFluxDG < FixedTotalFlux
    properties
    end
    
    methods
        function vT = evaluateOnDomain(prop, model, state)
            vT = sum(state.flux, 2);
            vT = vT(state.faces);
        end
    end
end