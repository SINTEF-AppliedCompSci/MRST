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
            for c = 1:ncomp
                mob = model.Components{c}.getComponentMobility(model, state);
                for ph = 1:nph
                    if ~isempty(mob{ph})
                        v{c, ph} = prop.faceUpstream(state, flag{ph}, -mob{ph}).*kgrad{ph};
                    end
                end
            end
        end
    end
end