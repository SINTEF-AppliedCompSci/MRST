classdef ComponentPhaseFlux < StateFunction
    properties

    end
    
    methods
        function cf = ComponentPhaseFlux(backend, upwinding)
            cf@StateFunction(backend);
            cf = cf.dependsOn({'PermeabilityPotentialGradient', 'FaceComponentMobility'});
        end

        function v = evaluateOnDomain(prop, model, state)
            ncomp = model.getNumberOfComponents;
            nph = model.getNumberOfPhases;
            [kgrad, compMob] = prop.getEvaluatedDependencies(state,...
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