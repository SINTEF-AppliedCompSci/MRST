classdef FixedTotalFluxDG < FixedTotalFlux
    properties
    end
    
    methods
        function vT = evaluateOnDomain(prop, model, state)
            vT = sum(state.flux, 2);
            [~, ~, ~, f] = model.disc.getCubature(find(model.operators.internalConn), 'face');
            vT = vT(f);
        end
    end
end