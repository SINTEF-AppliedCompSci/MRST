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
            if isfield(ds.data, s(1).subs)
                h = builtin('subsref', ds.data, s);
            else
                h = builtin('subsref', ds, s);
            end
        end
        
        function u = struct(v)
            u = v.data;
        end
        
        function u = value(v)
            u = value(v.data);
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
        
        function ok = isfield(ds, f)
            ok = isfield(ds.data, f);
        end
        
        function ok = structPropEvaluated(s, name)
            ok = ~isempty(s.data.(name));
        end
    end
end

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
