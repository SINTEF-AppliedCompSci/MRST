classdef EffectiveMixturePolymerViscMultiplier < StateFunction
    properties
    end

    methods
        function gp = EffectiveMixturePolymerViscMultiplier(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn({'polymer'}, 'state'); % check mechanism
        end

        function muWeffMult = evaluateOnDomain(prop, model, state)
            cp   = model.getProp(state, 'polymer');            
            fluid = model.fluid;
            mixpar = fluid.mixPar;
            cpbar   = cp/fluid.cmax;
            a = fluid.muWMult(fluid.cmax).^(1-mixpar);
            b = 1./(1 - cpbar + cpbar./a);
            % The viscosity multiplier only results from the polymer mixing.
            muWeffMult = b.*fluid.muWMult(cp).^mixpar;            
        end
    end
end