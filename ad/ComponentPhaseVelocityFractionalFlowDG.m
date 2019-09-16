classdef ComponentPhaseVelocityFractionalFlowDG < StateFunction
    properties
    end
    
    methods
        function cf = ComponentPhaseVelocityFractionalFlowDG(backend, upwinding)
            cf@StateFunction(backend);
            cf = cf.dependsOn({'TotalVelocity'});
%             cf = cf.dependsOn({'TotalVelocity', 'PhaseInterfacePressureDifferences', 'Transmissibility', 'Mobility', 'TotalMobility'});
            cf = cf.dependsOn({'ComponentMobility', 'Mobility', 'TotalMobility', 'GravityPermeabilityGradient'}, 'FlowPropertyFunctions');
        end

        function v = evaluateOnDomain(prop, model, state)
            ncomp = model.getNumberOfComponents;
            nph = model.getNumberOfPhases;
            
            
            vT = prop.getEvaluatedDependencies(state, 'TotalVelocity');
            
            compMob = model.getProp(state, 'ComponentMobility');
            gRhoKdz = model.getProp(state, 'GravityPermeabilityGradient');
            mob     = model.getProp(state, 'Mobility');
            
            mobT = 0;
            for i = 1:nph
                mobT = mobT + mob{i};
            end
            w = 1./mobT;
            kgrad = cell(1, nph);
            for i = 1:nph
                mobG = 0;
                for j = 1:nph
                    if i ~= j
                        mobG = mobG + mob{j}.*(gRhoKdz{i} - gRhoKdz{j});
                    end
                end
                kgrad{i} = w.*(vT + mobG);
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