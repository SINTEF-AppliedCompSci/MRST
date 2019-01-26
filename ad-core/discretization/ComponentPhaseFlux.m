classdef ComponentPhaseFlux < GridProperty & UpwindProperty
    properties

    end
    
    methods
        function cf = ComponentPhaseFlux(backend, upwinding)
            cf@GridProperty(backend);
            cf@UpwindProperty(upwinding)
        end

        
        function v = evaluateOnDomain(prop, model, state)
            ncomp = model.getNumberOfComponents;
            nph = model.getNumberOfPhases;
            [kgrad, flag] = model.getProps(state,...
                'PermeabilityPotentialGradient', 'PhaseUpwindFlag');
            v = cell(ncomp, nph);
            for i = 1:ncomp
                mob = model.Components{i}.getComponentMobility(model, state);
                for j = 1:nph
                    if ~isempty(mob{j})
                        v{i, j} = prop.faceUpstream(state, flag{j}, mob{i}).*kgrad{j};
                    end
                end
            end
        end
    end
end