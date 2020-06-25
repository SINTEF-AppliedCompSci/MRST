classdef MechPropertyFunctions < StateFunctionGrouping
    
    properties
        FaceNodeDisplacement
        Strain 
        Stress 
    end

    methods
        
        function props = MechPropertyFunctions(model)
            
            props@StateFunctionGrouping('MechProps');
            
            stress = Stress(model);
            strain = Strain(model);
            fndisp = FaceNodeDisplacement(model);
            
            props = props.setStateFunction('Stress', stress);
            props = props.setStateFunction('Strain', strain);
            props = props.setStateFunction('FaceNodeDisplacement', fndisp);
            
        end
    end
    
end
