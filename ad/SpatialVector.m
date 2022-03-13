classdef SpatialVector
    % Class for supporting ad as spatial vectors, for use in e.g.,
    % Discontinuous Galerkin discretizations
   
    properties
        vals
        dim
    end
    
    methods 
        function v = SpatialVector(varargin)
            if all(cellfun(@(v) isa(v, 'ADI'), varargin))
                v.dim       = nargin;
                v.vals      = cell(1,v.dim);
                [v.vals{:}] = varargin{:};
            else
                if all(cellfun(@(v) isa(v, 'double'), varargin))
                    sz = cellfun(@(v) size(v,2), varargin);
                    if all(sz == sz(1))
                        vals = varargin;
                        v.dim = nargin;
                        v.vals = vals;
                    else
                        vals = varargin{1}.*varargin{2};
                        v.dim = size(vals,2);
                        v.vals = mat2cell(vals, size(vals,1), ones(1,size(vals,2)));
                    end
                    return 
                end
                assert(nargin == 2)
                s = varargin{1};
                v = varargin{2};
                if ~isa(s, 'double')
                    v = SpatialVector(v,s);
                else
                    assert(any(strcmpi(class(v), getADtypes())));
                    s = mat2cell(s, size(s,1), ones(1,size(s,2)));
                    v = repmat({v}, 1, size(s,2));
                    v = cellfun(@times, s, v, 'uniformOutput', false);
                    v = SpatialVector(v{:});
                end
            end
        end
        
        function u = subsref(v, s)
            % Subscripted reference. Called for `h = u(v)`.
            if strcmp(s(1).type, '.')
                u = builtin('subsref',v,s);
            else
                switch s(1).type
                    case '()'
                        u = v;
                        subs = s(1).subs; 
                        switch numel(subs)
                            case 1
                                if ischar(subs{1}) && strcmpi(subs{1}, ':')
                                    subs{1} = 1:size(u,2);
                                end
                                u.vals = u.vals(subs{1});
                                u.dim  = numel(subs{1});
                            case 2
                                if ischar(subs{2}) && strcmpi(subs{2}, ':')
                                    subs{2} = 1:size(u,2);
                                end
                                u.vals = cellfun(@(u) u(subs{1}), u.vals(subs{2}), 'uniformOutput', false);
                                u.dim  = numel(subs{2});
                        end
                end
                if numel(s) > 1
                    % Recursively handle next operation
                    u = subsref(u, s(2:end));
                end
            end
        end
        
        function sz = size(v,dim)
            sz = [numel(value(v.vals{1})), numel(v.vals)];
            if nargin == 2
                sz = sz(dim);
            end
        end
        
        function u = uplus(v)
            u = v;
        end
        
        function u = uminus(v)
            u = v;
            u.vals = cellfun(@uminus, v.vals, 'uniformOutput', false);
        end
        
        function r = plus(u,v)
            r = u;
            if ~isa(u, 'SpatialVector')
                r = plus(v,u);
            else
                if ~isa(v, 'SpatialVector')
                    if isa(v, 'double')                    
                        sz = size(v);
                        v = mat2cell(v, sz(1), ones(1, sz(2)));
                        if numel(v) < numel(u.vals)
                            v = repmat(v, 1, size(u,2));
                        end
                    elseif isa(v, 'ADI')
                        v = repmat({v}, 1, u.dim);
                    end
                else
                    v = v.vals;
                end
                r.vals = cellfun(@plus, u.vals, v, 'uniformOutput', false);
            end
        end
        
        function r = minus(u,v)
            r = plus(u,-v);
        end
        
        function s = sum(u,dim)
            if nargin == 1
                dim = 2;
            end
            switch dim
                case 1
                    s = u;
                    s.vals = cellfun(@sum, u.vals, 'uniformOutput', false);
                case 2
                    s = 0;
                    for dNo = 1:u.dim
                        s = s + u.vals{dNo};
                    end
            end
            
        end
        
        function r = times(u,v)
            %assert(all(size(u) == size(v)));
            r = u;
            if ~isa(u, 'SpatialVector')
                r = times(v,u);
            else
                if ~isa(v, 'SpatialVector')
                    if isa(v, 'double')                    
                        sz = size(v);
                        v = mat2cell(v, sz(1), ones(1, sz(2)));
                        if numel(v) < numel(u.vals)
                            v = repmat(v, 1, size(u,2));
                        end
                    elseif isa(v, 'ADI')
                        v = repmat({v}, 1, u.dim);
                    end
                else
                    v = v.vals;
                end 
                r.vals = cellfun(@times, u.vals, v, 'uniformOutput', false);
            end
        end
        
        function r = dot(u,v)
            r = sum(times(u,v),2);
        end
        
        function r = cross(u,v)
            assert(u.dim == 3 && v.dim == 3)
            r = u;
            r.vals{1} =   u.vals{2}.*v.vals{3} - u.vals{3}.*v.vals{2} ;
            r.vals{2} = -(u.vals{1}.*v.vals{3} - u.vals{3}.*v.vals{1});
            r.vals{3} =   u.vals{1}.*v.vals{2} - u.vals{2}.*v.vals{1} ;
        end
        
    end
    
end

function ad = getADtypes()
    ad = {'ADI'};
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
