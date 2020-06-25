classdef BiotPropertyFunctions < StateFunctionGrouping
    
    properties
        Dilatation
        BasePoreVolume
    end

    methods
        function props = BiotPropertyFunctions(model)
            
            props@StateFunctionGrouping('BiotProps');
            
            dilatation = BiotCoupledDilatation(model);
            pv = PoreVolume(model);
            
            props = props.setStateFunction('Dilatation', dilatation);
            props = props.setStateFunction('BasePoreVolume', pv);
            
        end
    end
    
end
