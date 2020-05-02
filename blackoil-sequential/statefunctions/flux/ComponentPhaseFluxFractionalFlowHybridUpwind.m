classdef ComponentPhaseFluxFractionalFlowHybridUpwind < StateFunction
    % Hybrid upwind (for gravity + viscous flow)
    properties
    end
    
    methods
        function cf = ComponentPhaseFluxFractionalFlowHybridUpwind(backend, upwinding)
            cf@StateFunction(backend);
            cf = cf.dependsOn({'FaceComponentMobility', ...
                               'TotalFlux', ...
                               'PhaseInterfacePressureDifferences', ...
                               'Transmissibility', ...
                               'FaceMobility', ...
                               'FaceTotalMobility', ...
                               'FaceComponentMobilityGravity', ...
                               'FaceMobilityGravity', ...
                               'FaceTotalMobilityGravity' ...
                               });
        end

        function v = evaluateOnDomain(prop, model, state)
            ncomp = model.getNumberOfComponents;
            nph = model.getNumberOfPhases;
            [compMob, vT, G, T, mobT] = prop.getEvaluatedDependencies(state,...
                 'FaceComponentMobility', 'TotalFlux', 'PhaseInterfacePressureDifferences', 'Transmissibility', 'FaceTotalMobility');
            [compMob_g, fmob_g, mobT_g] = prop.getEvaluatedDependencies(state,...
                'FaceComponentMobilityGravity', 'FaceMobilityGravity', 'FaceTotalMobilityGravity');
            w = 1./mobT;
            w_g = 1./mobT_g;
            kgrad = cell(1, nph);
            kgrad_g = cell(1, nph);
            for i = 1:nph
                mobG = 0;
                for j = 1:nph
                    if i ~= j
                        mobG = mobG + fmob_g{j}.*(G{i} - G{j});
                    end
                end
                kgrad{i} = w.*vT;
                kgrad_g{i} = w_g.*T.*mobG;
            end
            v = cell(ncomp, nph);
            for c = 1:ncomp
                for ph = 1:nph
                    mob = compMob{c, ph};
                    mob_g = compMob_g{c, ph};
                    if ~isempty(mob)
                        v{c, ph} = mob.*kgrad{ph} + mob_g.*kgrad_g{ph};
                    end
                end
            end
        end
    end
end