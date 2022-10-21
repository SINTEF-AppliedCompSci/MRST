classdef WellboreBottomholeTemperature < StateFunction
    
    methods
        %-----------------------------------------------------------------%
        function cf = WellboreBottomholeTemperature(model)

            cf@StateFunction(model);
            cf = cf.dependsOn({'T'}, 'state');
            cf.label = 'T_{bh}';

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function bht = evaluateOnDomain(prop, model, state)
            
            T   = model.getProps(state, 'T');
            iic = model.getInletSegments();
            bht = T(iic);
            
        end
        %-----------------------------------------------------------------%
        
    end
    
end