classdef TransmissibilityDG < Transmissibility
    properties
        
    end
    
    methods
        function pp = TransmissibilityDG(model)
            pp@Transmissibility(model);
        end
        
        function T = evaluateOnDomain(prop, model, state)
            T = model.operators.T_all;
            [~, ~, ~, f] = model.disc.getCubature(find(model.operators.internalConn), 'face');
            T = T(f);
            if isfield(model.fluid, 'transMult')
                p = model.getProps(state, 'pressure');
                T = model.fluid.transMult(p).*T;
            end
        end
    end
end