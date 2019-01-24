classdef Transmissibility < GridProperty
    properties
        
    end
    
    methods
        function pp = Transmissibility(backend)
            pp@GridProperty(backend);
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