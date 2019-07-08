classdef ComponentPhaseFluxFractionalFlow < StateFunction
    properties
    end
    
    methods
        function cf = ComponentPhaseFluxFractionalFlow(backend, upwinding)
            cf@StateFunction(backend);
            cf = cf.dependsOn({'FaceComponentMobility', 'TotalFlux', 'PhaseInterfacePressureDifferences', 'Transmissibility', 'FaceMobility', 'FaceTotalMobility'});
        end

        function v = evaluateOnDomain(prop, model, state)
            ncomp = model.getNumberOfComponents;
            nph = model.getNumberOfPhases;
            [compMob, vT, G, T, fmob, mobT] = prop.getEvaluatedDependencies(state,...
                 'FaceComponentMobility', 'TotalFlux', 'PhaseInterfacePressureDifferences', 'Transmissibility', 'FaceMobility', 'FaceTotalMobility');
            w = 1./mobT;
            kgrad = cell(1, nph);
            for i = 1:nph
                mobG = 0;
                for j = 1:nph
                    if i ~= j
                        mobG = mobG + fmob{j}.*(G{i} - G{j});
                    end
                end
                kgrad{i} = w.*(vT + T.*mobG);
            end
            v = cell(ncomp, nph);
            for c = 1:ncomp
                for ph = 1:nph
                    mob = compMob{c, ph};
                    if ~isempty(mob)
                        v{c, ph} = mob.*kgrad{ph};
                    end
                end
            end
        end
    end
end