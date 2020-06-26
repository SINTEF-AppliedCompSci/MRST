classdef ConsistentDiv < StateFunction
    
    methods
        function gp = ConsistentDiv(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn('FaceNodeDisplacement');
            gp = gp.dependsOn('displacement', 'state');
        end
        
        function cdiv = evaluateOnDomain(prop, model, state)
            cdivop = model.operators.cdivop;
            [unf, uc] = model.getProps(state, 'FaceNodeDisplacement', 'displacement');
            cdiv = cdivop(unf, uc);
        end
    end
    
end
