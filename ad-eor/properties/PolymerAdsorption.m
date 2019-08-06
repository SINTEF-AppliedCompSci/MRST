classdef PolymerAdsorption < GridProperty
    properties
    end

    methods
        function gp = PolymerAdsorption(model, varargin)
            gp@GridProperty(model, varargin{:});
            gp = gp.dependsOn({'polymer', 'polymermax'}, 'state'); % check mechanism
        end

        function ads = evaluateOnDomain(prop, model, state)
            [cp, cpmax] = model.getProps(state, 'polymer', 'polymermax');
            fluid = model.fluid;
            ads  = effads(cp, cpmax, fluid);
        end
    end
end

function y = effads(cp, cpmax, f)
   if f.adsInx == 2
      y = f.ads(max(cp, cpmax));
   else
      y = f.ads(cp);
   end
end