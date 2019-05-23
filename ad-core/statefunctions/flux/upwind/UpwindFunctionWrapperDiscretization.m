classdef UpwindFunctionWrapperDiscretization < UpwindDiscretization
    % Simple wrapper
    properties (Access = protected)
        function_handle
    end
    
    methods
        function ufn = UpwindFunctionWrapperDiscretization(model)
            ufn@UpwindDiscretization(model);
            if isempty(model.operators)
                N = getNeighbourship(model.G);
                nc = model.G.cells.num;
                nf = size(N, 1);
                up = @(flag, x)faceUpstr(flag, x, N, [nf, nc]);
            else
                up = model.operators.faceUpstr;
            end
            assert(isa(up, 'function_handle'));
            ufn.function_handle = up;
        end
        
        function v = faceUpstream(wrapper, model, state, flag, cellvalue)
            v = wrapper.function_handle(flag, cellvalue);
        end
    end
end