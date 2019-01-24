classdef ComponentFlux < GridProperty
    properties

    end
    
    methods
        function v = evaluateOnDomain(prop, model, state)
            ncomp = numel(model.Components);
            [kgrad, mobf] = model.getProps(state, 'PermeabilityPotentialGradient');
            v = cell(1, ncomp);
            for i = 1:ncomp
                mob = model.Components{i}.getComponentMobility(model, state);
                v{i} = mobf.*kgrad;
            end
        end
    end
end