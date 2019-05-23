classdef PhaseFlux < StateFunction
    properties

    end
    
    methods
        function fm = PhaseFlux(model)
            fm@StateFunction(model);
            fm = fm.dependsOn({'FaceMobility', 'PermeabilityPotentialGradient'});
        end

        
        function v = evaluateOnDomain(prop, model, state)
            [mob, kgrad] = prop.getEvaluatedDependencies(state,...
                'FaceMobility', 'PermeabilityPotentialGradient');
            nph = numel(mob);
            v = cell(1, nph);
            for i = 1:nph
                v{i} = -mob{i}.*kgrad{i};
            end
        end
    end
end