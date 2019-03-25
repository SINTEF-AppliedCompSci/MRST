classdef ComponentPhaseFlux < GridProperty
    properties

    end
    
    methods
        function cf = ComponentPhaseFlux(backend, upwinding)
            cf@GridProperty(backend);
        end

        
        function v = evaluateOnDomain(prop, model, state)
            ncomp = model.getNumberOfComponents;
            nph = model.getNumberOfPhases;
            [kgrad, compMob] = model.getProps(state,...
                'PermeabilityPotentialGradient', 'FaceComponentMobility');
            v = cell(ncomp, nph);
            for c = 1:ncomp
                for ph = 1:nph
                    mob = compMob{c, ph};
                    if ~isempty(mob)
                        v{c, ph} = -mob.*kgrad{ph};
                    end
                end
            end
        end
    end
end