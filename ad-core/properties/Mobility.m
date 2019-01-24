classdef Mobility < GridProperty
    properties
    end
    
    methods
        function mob = evaluateOnDomain(prop, model, state)
            [mu, kr] = model.getProps(state, 'Viscosity', 'RelativePermeability');
            mob = cellfun(@(x, y) x./y, kr, mu, 'unif', false);
        end
    end
end