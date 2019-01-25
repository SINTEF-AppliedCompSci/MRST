classdef ComponentFlux < GridProperty & UpwindProperty
    properties

    end
    
    methods
        function cf = ComponentFlux(backend, upwinding)
            cf@GridProperty(backend);
            cf@UpwindProperty(upwinding)
        end

        
        function v = evaluateOnDomain(prop, model, state)
            ncomp = numel(model.Components);
            [kgrad, flag] = model.getProps(state,...
                'PermeabilityPotentialGradient', 'PhaseUpwindFlag');
            v = cell(1, ncomp);
            for i = 1:ncomp
                mob = model.Components{i}.getComponentMobility(model, state);
                v{i} = prop.faceUpstream(state, flag{i}, mob).*kgrad;
            end
        end
    end
end