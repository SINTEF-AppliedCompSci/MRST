classdef FluidBiotPropertyFunctions < StateFunctionGrouping
    
    properties 
        
        PoreVolume
        MassConsEquations
        
    end

    methods
        
        function props = MechBiotPropertyFunctions(model)
            
            props@StateFunctionGrouping('FluidBiotProps');
            
            pv   = PoreVolume(model);
            masseqs   = MassConsEquations(model);
            
            props = props.setStateFunction('MassConsEquations', masseqs);
            props = props.setStateFunction('PoreVolume', pv);
            
        end
    end
    
end
