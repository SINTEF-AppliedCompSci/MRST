classdef Mobility < GridProperty
    properties
    end
    
    methods
        function mob = evaluateOnGrid(prop, model, state)
            props = model.FlowPropertyFunctions;
            mu = props.getProperty(model, state, 'Viscosity');
            kr = props.getProperty(model, state, 'RelativePermeability');
            mob = cellfun(@(x, y) x./y, kr, mu, 'unif', false);
        end
    end
end