classdef SurfactantViscMultiplier < StateFunction
    properties
    end

    methods
        function gp = SurfactantViscMultiplier(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn('pressure', 'state');
        end

        function muWMults = evaluateOnDomain(prop, model, state)
            p = model.getProps(state, 'pressure');
            fluid = model.fluid;
            muWMults = fluid.muW(p)/fluid.muWr;
        end
    end
end