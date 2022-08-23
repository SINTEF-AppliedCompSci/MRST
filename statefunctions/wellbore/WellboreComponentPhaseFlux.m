classdef WellboreComponentPhaseFlux < StateFunction
    
    methods
        %-----------------------------------------------------------------%
        function cf = WellboreComponentPhaseFlux(model)

            cf@StateFunction(model);
            cf = cf.dependsOn({'massFlux'}, 'state');
            cf.label = 'V_{i,\alpha}^w';

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function v = evaluateOnDomain(prop, model, state)
            
            v = {state.massFlux};
            
        end
        %-----------------------------------------------------------------%
        
    end
    
end