classdef WellboreEffect < StateFunction
    
    methods
        %-----------------------------------------------------------------%
        function cf = WellboreEffect(model)

            cf@StateFunction(model);
            cf = cf.dependsOn({'HeatFlux'}, 'FlowPropertyFunctions');
            cf.label = 'q_{h,s}';

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function qh = evaluateOnDomain(prop, model, state)
            
            qh = model.getProps(state, 'HeatFlux');
            [~, iif] = model.getInletSegments();
            qh = model.parentModel.operators.fluxSum(qh(iif));
            
        end
        %-----------------------------------------------------------------%
        
    end
    
end