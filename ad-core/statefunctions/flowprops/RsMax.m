classdef RsMax < StateFunction
    properties
    end
    
    methods
        function gp = RsMax(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn('pressure', 'state');
        end
        
        function rsSat = evaluateOnDomain(prop, model, state)
            p = model.getProp(state, 'pressure');
            if model.disgas
                f = model.fluid;
                rsSat = prop.evaluateFunctionOnDomainWithArguments(f.rsSat, p);
            else
                rsSat = 0*p;
            end
        end
    end
end