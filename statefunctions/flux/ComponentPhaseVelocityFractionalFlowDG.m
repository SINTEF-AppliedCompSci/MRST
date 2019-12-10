classdef ComponentPhaseVelocityFractionalFlowDG < StateFunction
    properties
    end
    
    methods
        function cf = ComponentPhaseVelocityFractionalFlowDG(backend)
            cf@StateFunction(backend);
            cf = cf.dependsOn({'TotalVelocity'});
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
            for ph = 1:nph
                mobT = mobT + mob{ph};
            end
            w = 1./mobT;
            kgrad = cell(1, nph);
            for ph = 1:nph
                mobG = 0;
                for j = 1:nph
                    if ph ~= j
                        mobG = mobG + mob{j}.*(gRhoKdz{ph} - gRhoKdz{j});
                    end
                end
                kgrad{ph} = w.*(vT + mobG);
            end
            v = cell(ncomp, nph);
            for c = 1:ncomp
                for ph = 1:nph
                    cmob = compMob{c, ph};
                    if ~isempty(cmob)
                        v{c, ph} = cmob.*kgrad{ph};
                    end
                end
            end
        end
    end
end