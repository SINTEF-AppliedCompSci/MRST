classdef DynamicStruct < handle
    properties (Access = private)
        data = struct();
    end
    
    methods
        function ds = DynamicStruct(s)
            if nargin > 0
                ds.data = s;
            end
        end
        
        function h = subsref(ds,s)
            h = builtin('subsref', ds.data, s);
        end
        
        function u = subsasgn(u,varargin)
            u.data = builtin('subsasgn',u.data, varargin{:});
        end
        
        function disp(ds)
            fprintf('Dynamic struct with fields:\n');
            disp(ds.data);
        end
    end
end