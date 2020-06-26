classdef BiotFaceNodeDisplacement < StateFunction
    
    methods
        function gp = BiotFaceNodeDisplacement(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'displacement', 'lambdamech', 'biotpressure'}, 'state');
        end
        
        function fndisp = evaluateOnDomain(prop, model, state)
            
            fndispop = model.operators.facenodedispop;
            [u, p, lm] = model.getProps(state, 'displacement', 'biotpressure', 'lambdamech');
            fndisp = fndispop(u, p, lm);
            
        end
    end
end
