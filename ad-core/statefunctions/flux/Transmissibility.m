classdef Transmissibility < StateFunction
    properties
        
    end
    
    methods
        function pp = Transmissibility(model)
            pp@StateFunction(model);
            if isfield(model.fluid, 'transMult')
                pp = pp.dependsOn('pressure', 'state');
            end
        end
        
        function v = evaluateOnDomain(prop, model, state)
            v = model.operators.T;
            if isfield(model.fluid, 'transMult')
                p = model.getProps(state, 'pressure');
                v = model.fluid.transMult(p).*T;
            end
        end
    end
end