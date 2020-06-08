classdef MechBiotPropertyFunctions < StateFunctionGrouping
    
    properties
        Strain 
        Stress 
        BiotGradP
        Dilation
        MomentumEquations
    end

    methods
        function props = MechBiotPropertyFunctions(model)
            
            props@StateFunctionGrouping('MechBiotProps');
            
            biotgradp   = BiotGradP(model);
            dilation    = BiotCoupledDilatation(model);
            stress      = Stress(model);
            strain      = Strain(model);
            momentumeqs = MomentumEquations(model);
            
            props = props.setStateFunction('Stress', stress);
            props = props.setStateFunction('Strain', strain);
            props = props.setStateFunction('BiotGradP', biotgradp);
            props = props.setStateFunction('Dilation', dilation);
            props = props.setStateFunction('MomentumEquations', momentumeqs);
            
        end
    end
    
end
