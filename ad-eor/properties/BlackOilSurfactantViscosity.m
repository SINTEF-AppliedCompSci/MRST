classdef BlackOilSurfactantViscosity < BlackOilViscosity
    properties
    end

    methods
        function gp = BlackOilSurfactantViscosity(model, varargin)
            gp@BlackOilViscosity(model, varargin{:});
            gp = gp.dependsOn('pressure', 'state');
        end

        function mu = evaluateOnDomain(prop, model, state)
           cs = model.getProps(state, 'surfactant');
           p = model.getProps(state, 'pressure');
           mu = prop.evaluateOnDomain@BlackOilViscosity(model, state);
           mu{1} = model.fluid.muWSft(cs);
           muWMults = model.fluid.muW(p)/model.fluid.muWr;
           
           mu{1} = mu{1}.*muWMults;
        end
    end
end