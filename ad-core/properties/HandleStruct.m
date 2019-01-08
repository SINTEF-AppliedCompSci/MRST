classdef HandleStruct < handle
    properties (Access = private)
        data = struct();
    end
    
    methods
        function ds = HandleStruct(s)
            if nargin > 0
                ds.data = s;
            end
        end
        
        function s = reduce(ds)
            s = ds.data;
            flds = fieldnames(s);
            for i = 1:numel(flds)
                f = flds{i};
                if isa(s.(f), 'ADI')
                    s.(f) = value(s.(f));
                elseif iscell(s.(f))
                    s.(f) = cellfun(@value, s.(f), 'UniformOutput', false);
                    s.(f) = [s.(f){:}];
                end
            end
        end
        
        function h = subsref(ds,s)
            if strcmp(s(1).type, '.') && ...
               ischar(s(1).subs) && ...
               isfield(ds.data, s(1).subs)
                h = builtin('subsref', ds.data, s);
            else
                h = builtin('subsref', ds, s);
            end
        end
        
        function u = struct(v)
            u = v.data;
        end
        
        function u = subsasgn(u,varargin)
            u.data = builtin('subsasgn',u.data, varargin{:});
        end
        
        function disp(ds)
            fprintf('HandleStruct with fields:\n');
            disp(ds.data);
        end
        
        function ds2 = copy(ds)
            % Deep copy
            ds2 = HandleStruct(ds.data);
        end
    end
end