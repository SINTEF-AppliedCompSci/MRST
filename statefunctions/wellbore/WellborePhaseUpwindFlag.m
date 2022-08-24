classdef WellborePhaseUpwindFlag < StateFunction
    
    methods
        %-----------------------------------------------------------------%
        function cf = WellborePhaseUpwindFlag(model)

            cf@StateFunction(model);
            cf = cf.dependsOn('ComponentPhaseFlux');

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function flag = evaluateOnDomain(prop, model, state)
            
            v = prop.getEvaluatedDependencies(state,...
                        'ComponentPhaseFlux');
            flag = cellfun(@(v) value(v) >= 0, v, 'UniformOutput', false);
            
        end
        %-----------------------------------------------------------------%
        
    end
    
end