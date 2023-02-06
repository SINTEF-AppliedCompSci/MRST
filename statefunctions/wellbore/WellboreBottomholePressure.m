classdef WellboreBottomholePressure < StateFunction
    
    methods
        %-----------------------------------------------------------------%
        function cf = WellboreBottomholePressure(model)

            cf@StateFunction(model);
            cf = cf.dependsOn({'pressure'}, 'state');
            cf.label = 'p_{bh}';

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function bhp = evaluateOnDomain(prop, model, state)
            
            p   = model.getProps(state, 'pressure');
            iic = model.getInletSegments();
            bhp = p(iic);
            
        end
        %-----------------------------------------------------------------%
        
    end
    
end