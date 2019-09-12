classdef PolymerViscMultiplier < StateFunction
    properties
    end

    methods
        function gp = PolymerViscMultiplier(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PolymerAdsorption'});
            gp = gp.dependsOn({'EffectiveMixturePolymerViscMultiplier'});
        end

        function muWMultp = evaluateOnDomain(prop, model, state)
            ads = model.getProp(state, 'PolymerAdsorption');
            muWeffMult = model.getProp(state, 'EffectiveMixturePolymerViscMultiplier');
            fluid = model.fluid;
            permRed = 1 + ((fluid.rrf - 1)./fluid.adsMax).*ads;
            muWMultp  = muWeffMult.*permRed;
        end
    end
end