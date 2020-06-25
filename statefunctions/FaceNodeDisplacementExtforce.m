classdef FaceNodeDisplacementExtforce < StateFunction
    
    methods
        function gp = FaceNodeDisplacementExtforce(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'displacement', 'lambdamech', 'pressure', 'extforce'}, 'state');
        end
        
        function fndisp = evaluateOnDomain(prop, model, state)
            
            fndispop = model.operators.facenodedispop;
            [u, p, lm, ef] = model.getProps(state, 'displacement', 'pressure', 'lambdamech', 'extforce');
            fndisp = fndispop(u, p, lm, ef);
            
        end
    end
end
