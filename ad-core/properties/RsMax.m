classdef RsMax < GridProperty
    properties
    end
    
    methods
        function rsSat = evaluateOnGrid(prop, model, state)
            p = model.getProp(state, 'pressure');
            if model.disgas
                f = model.fluid;
                rsSat = prop.evaluateFunctionOnGrid(f.rsSat, p);
            else
                rsSat = 0*p;
            end
        end
    end
end