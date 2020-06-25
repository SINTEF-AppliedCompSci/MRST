classdef BiotPropertyFunctions < StateFunctionGrouping
    
    properties
        Dilation
    end

    methods
        function props = BiotPropertyFunctions(model)
            
            props@StateFunctionGrouping('BiotProps');
            
            dilation = BiotCoupledDilatation(model);
            
            props = props.setStateFunction('Dilation', dilation);
            
        end
    end
    
end
