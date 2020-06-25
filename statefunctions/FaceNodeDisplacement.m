classdef FaceNodeDisplacement < StateFunction
    
    methods
        function gp = FaceNodeDisplacement(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'displacement', 'lambdamech'}, 'state');
        end
        
        function fndisp = evaluateOnDomain(prop, model, state)
            
            fndispop = model.operators.facenodedispop;
            [u, lm] = model.getProps(state, 'displacement', 'lambdamech');
            fndisp = fndispop(u, lm);
            
        end
    end
end
