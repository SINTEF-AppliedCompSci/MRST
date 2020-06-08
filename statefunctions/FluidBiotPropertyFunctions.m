classdef FluidBiotPropertyFunctions < StateFunctionGrouping
    
    properties 
        
        PoreVolume
        Density
        
    end

    methods
        
        function props = FluidBiotPropertyFunctions(model)
            
            props@StateFunctionGrouping('FluidBiotProps');
            
            pv   = PoreVolume(model);
            density   = PoreVolume(model);
            props = props.setStateFunction('PoreVolume', pv);
            props = props.setStateFunction('Density', density);
            
        end
    end
    
end
