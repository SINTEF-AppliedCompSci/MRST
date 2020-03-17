classdef DiagonalSubset < DiagonalJacobian
    % Structured subset of a diagonal jacobian 
    properties
        map % Map to the underlying DiagonalJacobian representation. Two DiagonalSubsets of the same map can be multiplied together, etc.
        parentSubset
    end
    
    methods
        function D = DiagonalSubset(d, dims, map, subset, parentSubset)
            if nargin == 0
                return
            end
            D.diagonal = d;
            D.dim = dims;
            D.map = map;
            if nargin > 3
                if islogical(subset)
                    subset = find(subset);
                end
                D.subset = subset;
                if nargin > 4
                    D.parentSubset = parentSubset;
                end
            end
        end

        function u = subsasgn(u, s, v)
            if strcmp(s(1).type, '.')
                u = builtin('subsasgn',u,s,v);
            else
                assert(false, 'Subsasgn not supported for diagonal subset.');
            end
        end
        
        function [x, D] = diagMult(v, x, D)
            if any(x.diagonal(:))
                x.diagonal = bsxfun(@times, x.diagonal, v);
            else
                x = x.toZero();
            end
        end
        
        function [x, D1, D2] = diagProductMult(v1, v2, x, y, D1, D2)
            persistent allow_implicit;
            if isempty(allow_implicit)
                allow_implicit = ~verLessThan('matlab','9.1');
            end
            numx = isnumeric(x);
            if numx || isnumeric(y)
                [x, D2] = diagMult(v2, x, D2);
                [y, D1] = diagMult(v1, y, D1);
                x = x + y;
            else
                if isempty(x.diagonal)
                    [x, D1] = diagMult(v1, y, D1);
                elseif isempty(y.diagonal)
                    [x, D2] = diagMult(v2, x, D2);
                else
                    if allow_implicit
                        x.diagonal = x.diagonal.*v2 + y.diagonal.*v1;
                    else
                        x.diagonal = bsxfun(@times, x.diagonal, v2) + ...
                                     bsxfun(@times, y.diagonal, v1);
                    end
                end
            end
        end

        function out = matrixDims(D, n)
            if isempty(D.subset)
                ni = size(D.map, 1);
            else
                ni = size(D.subset, 1)/size(D.map, 2);
            end
            out = [ni, prod(D.dim)];
            
            if nargin == 1
                return
            end
            
            out = out(n);
        end

        function isEqual = subsetsEqual(x, y)
            xSub = isa(x, 'DiagonalSubset');
            ySub = isa(y, 'DiagonalSubset');
            if xSub && ySub
                if all(size(x.map) == size(y.map)) && all(all(x.map == y.map))
                    isEqual = subsetsEqual@DiagonalJacobian(x, y);
                else
                    isEqual = false;
                end
            else
                isEqual = false;
            end
        end
        
        function [I, J, V, imax, jmax] = getSparseBlocks(D)
            n = size(D.diagonal, 1);
            m = D.dim(2);
            nmap = size(D.map, 2);
            imax = n;
            jmax = prod(D.dim);
            if isempty(D.subset)
                nval = size(D.diagonal, 1);
            else
                nval = numel(D.subset)/nmap;
            end
            if size(D.map, 1) == 1
                D.map = repmat(D.map, nval, 1);
            end
            I = repmat((1:n)', 1, nmap*m);
            jmap = D.map;
            if ~isempty(D.subset)
                jmap = jmap(D.subset, :);
            end
            if ~isempty(D.parentSubset)
                jmap = D.parentSubset(jmap);
            end
            V = D.diagonal;
            J = jmap(:, rldecode((1:nmap)', m));
            tmp = (0:m-1)*D.dim(1);
            J = bsxfun(@plus, J, repmat(tmp, 1, nmap));
            V(all(jmap == 0, 2), :) = 0;
        end
        
        function u = repmat(u, varargin)
            u = repmat(u.sparse(), varargin{:});
        end
    end
end

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
