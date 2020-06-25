classdef MechPropertyFunctions < StateFunctionGrouping
    
    properties
        Strain 
        Stress 
    end

    methods
        
        function props = MechPropertyFunctions(model)
            
            props@StateFunctionGrouping('MechProps');
            
            stress       = Stress(model);
            strain       = Strain(model);
            
            props = props.setStateFunction('Stress', stress);
            props = props.setStateFunction('Strain', strain);
            
        end
    end
    
end
