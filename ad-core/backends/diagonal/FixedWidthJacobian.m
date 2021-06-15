classdef FixedWidthJacobian < DiagonalJacobian
    % Structured subset of a diagonal jacobian 
    properties
        map % Map to the underlying DiagonalJacobian representation. Two FixedWidthJacobians of the same map can be multiplied together, etc.
        mapName = '';
        parentSubset
    end
    
    methods
        function D = FixedWidthJacobian(d, dims, map, subset, parentSubset, useMex, useRowMajorMemory, mapName)
            if nargin == 0
                return
            end
            D.diagonal = d;
            D.dim      = dims;
            D.map      = map;
            D.useMex   = D.useMex;
            D.rowMajor = D.rowMajor;
            if nargin > 3
                if islogical(subset)
                    subset = find(subset);
                end
                D.subset = subset;
                if nargin > 4
                    D.parentSubset = parentSubset;
                    if nargin > 5
                        D.useMex = useMex;
                        if nargin > 6
                            D.rowMajor = useRowMajorMemory;
                            if nargin > 7
                                D.mapName = mapName;
                            end
                        end
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
            xSub = isa(x, 'FixedWidthJacobian');
            ySub = isa(y, 'FixedWidthJacobian');
            if xSub && ySub
                if x.mapsEqual(y)
                    isEqual = subsetsEqual@DiagonalJacobian(x, y);
                else
                    isEqual = false;
                end
            else
                isEqual = false;
            end
        end
        
        
        function isEqual = subsetsEqualNoZeroCheck(x, y)
            xSub = isa(x, 'FixedWidthJacobian');
            ySub = isa(y, 'FixedWidthJacobian');
            if xSub && ySub
                if x.mapsEqual(y)
                    isEqual = subsetsEqualNoZeroCheck@DiagonalJacobian(x, y);
                else
                    isEqual = false;
                end
            else
                isEqual = false;
            end
        end
        
        function isEqual = mapsEqual(x, y)
            if isempty(x.mapName)
                isEqual = all(size(x.map) == size(y.map)) && all(all(x.map == y.map));
            else
                isEqual = strcmp(x.mapName, y.mapName);
            end
        end
        
        function s = sparse(D)
            s = sparse@DiagonalJacobian(D);
        end
        
        function [I, J, V, imax, jmax] = getSparseBlocks(D)
            [I, J, V, imax, jmax] = getSparseBlocksFixedWidth(D);
        end
        
        function u = repmat(u, varargin)
            u = repmat(u.sparse(), varargin{:});
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
