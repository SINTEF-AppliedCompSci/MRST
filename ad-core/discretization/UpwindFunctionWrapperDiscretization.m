classdef UpwindFunctionWrapperDiscretization < UpwindDiscretization
    % Simple wrapper
    properties (Access = protected)
        function_handle
    end
    
    methods
        function ufn = UpwindFunctionWrapperDiscretization(fn)
            ufn@UpwindDiscretization();
            assert(isa(fn, 'function_handle'));
            ufn.function_handle = fn;
        end
        
        function v = faceUpstream(wrapper, state, flag, cellvalue)
            v = wrapper.function_handle(flag, cellvalue);
        end
    end
end