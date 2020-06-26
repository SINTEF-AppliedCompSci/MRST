classdef MechPropertyFunctions < StateFunctionGrouping
    
    properties
        FaceNodeDisplacement
        ConsistentDiv
        Strain 
        Stress 
    end

    methods
        
        function props = MechPropertyFunctions(model)
            
            props@StateFunctionGrouping('MechProps');
            
            stress = Stress(model);
            strain = Strain(model);
            fndisp = FaceNodeDisplacement(model);
            cdiv   = ConsistentDiv(model);
            
            props = props.setStateFunction('Stress', stress);
            props = props.setStateFunction('Strain', strain);
            props = props.setStateFunction('FaceNodeDisplacement', fndisp);
            props = props.setStateFunction('ConsistentDiv', cdiv);
            
        end
    end
    
end
