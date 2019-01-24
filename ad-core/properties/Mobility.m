classdef Mobility < GridProperty
    properties
    end
    
    methods
        function mob = evaluateOnDomain(prop, model, state)
            [mu, kr] = model.getProps(state, 'Viscosity', 'RelativePermeability');
            mob = cellfun(@(x, y) x./y, kr, mu, 'UniformOutput', false);
            if isfield(model.fluid, 'tranMultR')
                % Pressure dependent mobility multiplier 
                p = model.getProp(state, 'pressure');
                mult = model.fluid.tranMultR(p);
                mob = cellfun(@(x) x.*mult, mob, 'UniformOutput', false);
            end
        end
    end
end