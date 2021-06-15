classdef DiagonalJacobian
    % Diagonal representation of a Jacobian
    properties
        diagonal % Dense matrix of diagonal derivatives
        dim % Vector: First dimension is the number of variables in block, while the second is the number of columns
        subset % Indices corresponding to the subset (if empty, class contains the full set)
        useMex = false;
        rowMajor = false;
        allowImplicitExpansion = mrstPlatform('octave') || mrstPlatform('major') >= 9
    end

    methods
        function D = DiagonalJacobian(d, dim, subset, useMex, useRowMajorMemory)
            if nargin == 0
                return
            end
            if nargin < 3
                subset = [];
            elseif islogical(subset)
                subset = find(subset);
            end
            D.diagonal = d;
            D.dim = dim;
            D.subset = subset;
            if nargin > 3
                D.useMex = useMex;
                if nargin > 4
                    D.rowMajor = useRowMajorMemory;
                end
            end
        end
        
        function u = setDiagonal(u, subs, diagonals)
            if nargin < 2
                u.diagonal = diagonals;
            else
                if u.rowMajor
                    u.diagonal(:, subs) = diagonals;
                else
                    u.diagonal(subs, :) = diagonals;
                end
            end
        end
        
        function n = getNumberOfDiagonals(u)
            n = u.dim(2);
        end
        
        function d = getDiagonalByIndex(u, index)
            d = u.diagonal;
            if ~isempty(d)
                if u.rowMajor
                    d = d(index, :);
                else
                    d = d(:, index);
                end
            end
        end

        function D = getDiagonal(u, subs)
            if nargin == 1
                D = u.diagonal;
            else
                if u.rowMajor
                    D = u.diagonal(:, subs);
                else
                    D = u.diagonal(subs, :);
                end
            end
        end

        function sub = getSubset(u)
            d = 1 + u.rowMajor;
            if u.isZero
                sub = zeros(size(u.diagonal, d), 1);
            elseif isempty(u.subset)
                sub = reshape(1:size(u.diagonal, d), [], 1);
            else
                sub = u.subset;
            end
        end

        function u = toZero(u, n)
            if nargin == 1
                % The function is not taking a subset, just converting an
                % existing subset to zero type if present.
                n = size(u.diagonal, 1+u.rowMajor);
                if ~isempty(u.subset)
                    u.subset = zeros(n, 1);
                end
            else
                % We recieved a subset - set it to zeros to indicate that
                % it is arbitrary, but present.
                u.subset = zeros(n, 1);
            end
            if u.rowMajor
                u.diagonal = zeros(0, n);
            else
                u.diagonal = zeros(n, 0);
            end
        end

        function u = expandZero(u)
            if u.rowMajor
                u.diagonal = zeros(u.dim(2), size(u.diagonal, 2));
            else
                u.diagonal = zeros(size(u.diagonal, 1), u.dim(2));
            end
        end
        
        function isz = isZero(u)
            isz = size(u.diagonal, 2 - u.rowMajor) == 0;
        end
        
        function [I, J, V, imax, jmax] = getSparseBlocks(D, ioffset, joffset)
            rm = D.rowMajor;
            if rm
                [m, imax] = size(D.diagonal);
            else
                [imax, m] = size(D.diagonal);
            end
            jmax = prod(D.dim);
            if m == 0
                I = zeros(0, 1); J = zeros(0, 1); V = zeros(0, 1);
                return
            end
            if isempty(D.subset)
                J = 1:jmax;
                subs = 1:D.dim(1);
            else
                subs = D.subset;
                
                if m == 1
                    J = subs;
                else
                    J = bsxfun(@plus, repmat(subs, 1, m), (0:m-1)*D.dim(1));
                end
            end
            I = repmat(reshape(1:numel(subs), [], 1), m, 1);
            if rm
                V = D.diagonal';
            else
                V = D.diagonal;
            end
            V(subs == 0, :) = 0;
            if nargin > 1
                I = I + ioffset;
                if nargin > 2
                    J = J + joffset;
                end
            end
        end
        
        function [I, J, V, imax, jmax] = getSparseArguments(D, ioffset, joffset)
            [I, J, V, imax, jmax] = getSparseBlocks(D);
            act = V ~= 0;
            I = reshape(I(act), [], 1);
            J = reshape(J(act), [], 1);
            V = reshape(V(act), [], 1);
            if nargin > 1
                I = I + ioffset;
                if nargin > 2
                    J = J + joffset;
                end
            end
        end
        
        function s = sparse(D)
            % Convert diagonal representation into a sparse representation
            if D.useMex && isempty(D.subset) && ...
              ~isempty(D.diagonal) && ~isa(D, 'FixedWidthJacobian')
                % MEX acceleration only for the simplest case
                s = mexDiagonalSparse(D.diagonal, D.subset, D.dim, D.rowMajor);
            else
                [I, J, V, n, m] = D.getSparseArguments();
                s = sparse(I, J, V, n, m);
            end
        end
        
        function s = double(D)
            % Convert to sparse
            s = D.sparse();
        end
        function s = full(D)
            % Convert to full matrix - may be very expensive!
            s = full(D.sparse());
        end
        
        function isEqual = subsetsEqualNoZeroCheck(x, y, xsubset, ysubset)
            % Check if subsets are equal (ignoring zeros in subset as
            % "always fits")
            if nargin == 2
                xsubset = x.subset;
                ysubset = y.subset;
            end
            chx = isempty(xsubset);
            chy = isempty(ysubset);
            if chx && chy
                isEqual = true;
            elseif chx
                if numel(ysubset) == x.dim(1)
                    tmp = (1:x.dim(1))';
                    isEqual = all(y.subset(ysubset ~= 0) == tmp(ysubset ~= 0));
                else
                    isEqual = false;
                end
            elseif chy
                if numel(xsubset) == x.dim(1)
                    tmp = (1:x.dim(1))';
                    isEqual = all(xsubset(xsubset ~= 0) == tmp(xsubset ~= 0));
                else
                    isEqual = false;
                end
            else
                isEqual = numel(xsubset) == numel(ysubset) && ...
                    all(xsubset  == ysubset | xsubset == 0 | y.subset == 0);
            end
        end

        function isEqual = subsetsEqual(x, y)
            if x.isZero || y.isZero
                isEqual = true;
            else
                isEqual = subsetsEqualNoZeroCheck(x, y);
            end
        end
        
        function u = repmat(u, varargin)
            u = repmat(sparse(u), varargin{:});
        end

        function x = subsetPlus(x, v, subs)
            % Increment a subset of the Jacobian in place. This can be
            % faster than using subsasgn that may produce intermediate
            % class instances, for example
            % x(subs) = x(subs) + y(subs)
            % will end up as
            % tmp1 = x(subs);
            % tmp2 = y(subs);
            % tmp3 = tmp1 + tmp2;
            % x(subs) = tmp3;
            if isa(x, 'DiagonalJacobian')
                if isa(v, 'DiagonalJacobian')
                    if v.isZero
                        return
                    end
                    if x.isZero
                        x = x.sparse();
                        x(subs, :) = x(subs, :) + sparse(v);
                        return
                    end
                    % Get subset, check individual values
                    s = x.getSubset();
                    if subsetsEqualNoZeroCheck(x, v, s(subs), v.subset)
                        if x.rowMajor
                            x.diagonal(:, subs) = x.diagonal(:, subs) + v.diagonal;
                        else
                            x.diagonal(subs, :) = x.diagonal(subs, :) + v.diagonal;
                        end
                        % x.subset = DiagonalJacobian.treatSubset(x.subset, y.subset);
                    else
                        x = x.sparse();
                        v = v.sparse();
                        x(subs, :) = x(subs, :) + v;
                    end
                elseif issparse(v)
                    x = x.sparse();
                    x(subs, :) = x(subs, :) + v;
                end
            else
                x(subs, :) = x(subs, :) + v;
            end
        end
        function x = subsetMinus(x, subs, v)
            x = subsetPlus(x, -v, subs);
        end
        function u = subsref(u, s)
            if strcmp(s(1).type, '.')
                u = builtin('subsref',u,s);
            else
                switch s(1).type
                    case '()'
                        ix = s(1).subs{1};
                        if u.isZero
                            if ~ischar(ix)
                                u = u.toZero(numel(ix));
                            end
                        elseif numel(s(1).subs) == 2 && ischar(s(1).subs{2})
                            if any(ix)
                                if u.rowMajor
                                    u.diagonal = u.diagonal(:, ix);
                                else
                                    u.diagonal = u.diagonal(ix, :);
                                end
                            else
                                u = u.toZero(0);
                            end
                            
                            if isempty(u.subset)
                                if ischar(ix)
                                    % Do nothing
                                else
                                    u.subset = reshape(ix, [], 1);
                                end
                            else
                                subs = u.getSubset();
                                u.subset = subs(ix);
                            end
                        else
                            u = u.sparse();
                            u = builtin('subsref', u, s);
                        end
                        
                        if numel(s) > 1
                            % Recursively handle next operation
                            u = subsref(u, s(2:end));
                        end
                    case '{}'
                        error('Operation not supported');
                end
            end
        end

        function u = subsasgn(u, s, v)
            if strcmp(s(1).type, '.')
                u = builtin('subsasgn',u,s,v);
            else
                switch s(1).type
                    case '()'
                        uD = isa(u, 'DiagonalJacobian');
                        vD = isa(v, 'DiagonalJacobian');

                        if uD && vD
                            uz = u.isZero();
                            vz = v.isZero();
                            if uz && vz
                                % Both are zero - nothing to do
                                return;
                            elseif vz
                                % We are inserting a zero Jacobian
                                if u.rowMajor
                                    u.diagonal(:, s.subs{1}) = 0;
                                else
                                    u.diagonal(s.subs{1}, :) = 0;
                                end
                            elseif uz
                                % We are inserting into a zero Jacobian
                                vsub = v.getSubset();
                                
                                if isempty(u.subset)
                                    if islogical(s.subs{1})
                                        s.subs{1} = find(s.subs{1});
                                    end
                                    doZero = ~ischar(s.subs{1}) && ~all(s.subs{1} == vsub);
                                else
                                    doZero = true;
                                end
                                u = u.expandZero();
                                if doZero
                                    u.subset = u.getSubset();
                                    u.subset(s.subs{1}) = vsub;
                                end
                                if u.rowMajor
                                    u.diagonal(:, s.subs{1}) = v.diagonal;
                                else
                                    u.diagonal(s.subs{1}, :) = v.diagonal;
                                end
                            else
                                usub = u.subset;
                                vsub = v.subset;
                                if subsetsEqualNoZeroCheck(u, v, usub, vsub)
                                    % If the subsets match directly, we are
                                    % ready
                                    allowDiag = true;
                                elseif isempty(usub)
                                    % Check that we are inserting into
                                    % right position for diagonal to make
                                    % sense
                                    allowDiag = ~isempty(vsub) && u.compareIndices(u, s.subs{1}, vsub);
                                else
                                    if isempty(vsub)
                                        vsub = (1:v.dim(1))';
                                    end
                                    % u.subset is not zero
                                    allowDiag = u.compareIndices(u, s.subs{1}, usub(s.subs{1})) &&...
                                                u.compareIndices(u, s.subs{1}, vsub);
                                end

                                if allowDiag
                                    % We can safely insert diagonals. The
                                    % next step is to actually do it.
                                    if u.rowMajor
                                        u.diagonal(:, s.subs{1}) = v.diagonal;
                                    else
                                        u.diagonal(s.subs{1}, :) = v.diagonal;
                                    end
                                    if ~isempty(usub) && ~isempty(vsub)
                                        % Handle zero (wild card) subsets
                                        u.subset(s.subs{1}) = max(usub(s.subs{1}), vsub);
                                    end
                                else
                                    u = u.sparse();
                                    v = v.sparse();
                                    u(s.subs{1}, :) = v;
                                end
                            end
                        elseif uD
                            % v is double
                            if numel(v) == 1
                                if u.rowMajor
                                    u.diagonal(:, s.subs{1}) = v;
                                else
                                    u.diagonal(s.subs{1}, :) = v;
                                end
                            else
                                u = u.sparse();
                                % Re-call to insert sparse matrices
                                u = subsasgn(u, s, v);
                            end
                        elseif vD
                            % u is double
                            v = v.sparse();
                            % Re-call to insert sparse matrices
                            u = subsasgn(u, s, v);
                        end
                    case '{}'
                        error('Operation not supported');
                end
            end
        end
        
        
        function x = uminus(x)
            x.diagonal = -x.diagonal;
        end
        
        function x = minus(x, y)
            x = plus(x, -y);
        end
        
        function x = power(x, e)
            x.diagonal = x.diagonal.^e;
        end


        function u = rdivide(u,v)
            % Right matrix divide: `h=u/v`
            if isscalar(v)
                if ~u.isZero()
                    if u.rowMajor
                        u.diagonal = u.diagonal./v';
                    else
                        u.diagonal = u.diagonal./v;
                    end
                    
                end
            else
                u = u.sparse();
                if isa(v, 'DiagonalJacobian')
                    v = v.sparse();
                end
                u = u./v;
            end
        end
        function x = plus(x, y)
            if DiagonalJacobian.isAllZeros(x)
                x = y;
                return
            end
            if DiagonalJacobian.isAllZeros(y)
                return
            end
            xD = isa(x, 'DiagonalJacobian');
            if xD && x.isZero
                x = y;
                return
            end
            yD = isa(y, 'DiagonalJacobian');
            if yD && y.isZero 
                return
            end

            if xD && yD
                if subsetsEqualNoZeroCheck(x, y)
                    x.diagonal = x.diagonal + y.diagonal;
                    x.subset = DiagonalJacobian.treatSubset(x.subset, y.subset);
                else
                    x = x.sparse() + y.sparse();
                end
                return
            end
            
            if xD
                if ~DiagonalJacobian.isAllZeros(y)
                    x = x.sparse();
                    x = plus(x, y);
                end
            else
                x = plus(y, x);
            end
        end
        
        function x = times(x, y)
            if DiagonalJacobian.isAllZeros(x)
                return
            end
            
            if DiagonalJacobian.isAllZeros(y)
                x = y;
                return
            end
            
            xD = isa(x, 'DiagonalJacobian');
            yD = isa(y, 'DiagonalJacobian');
            if xD && yD
                if subsetsEqualNoZeroCheck(x, y)
                    x.diagonal = x.diagonal.*y.diagonal;
                    x.subset = DiagonalJacobian.treatSubset(x.subset, y.subset);
                else
                    x = x.sparse().*y.sparse();
                end
                return
            end
            
            if xD
                if isscalar(y)
                    x.diagonal = x.diagonal*y;
                else
                    x = x.sparse();
                    x = x.*y;
                end
            else
                x = times(y, x);
            end
        end
        
        function x = mtimes(x, y)
            if DiagonalJacobian.isAllZeros(x)
                if isscalar(x)
                    % Scalar multiply by zero - still diagonal
                    y = y.toZero();
                    x = y;
                else
                    % Assume valid matrix
                    dx = matrixDims(x);
                    dy = matrixDims(y);
                    x = sparse([], [], [], dx(1), dy(2));
                end
            elseif DiagonalJacobian.isAllZeros(y)
                if isscalar(y)
                    x = x.toZero();
                else
                    dx = matrixDims(x);
                    dy = matrixDims(y);
                    x = sparse([], [], [], dx(1), dy(2));
                end
            else
                xD = isa(x, 'DiagonalJacobian');
                yD = isa(y, 'DiagonalJacobian');
                if xD && yD
                    if subsetsEqualNoZeroCheck(x, y)
                        x.diagonal = x.diagonal.*y.diagonal;
                        x.subset = DiagonalJacobian.treatSubset(x.subset, y.subset);
                    else
                        x = x.sparse().*y.sparse();
                    end
                elseif xD
                    if isscalar(y)
                        x.diagonal = x.diagonal*y;
                    else
                        x = x.sparse();
                        x = x*y;
                    end
                else
                    if isscalar(x)
                        y.diagonal = x*y.diagonal;
                        x = y;
                    else
                        y = y.sparse();
                        x = x*y;
                    end
                end
            end
        end
        
        function varargout = matrixDims(D, n)
            dims = [size(D.diagonal, 1+D.rowMajor), prod(D.dim)];
            if nargout == 1
                varargout{1} = dims;
                if nargin > 1
                    varargout{1} = varargout{1}(n);
                end
            else
                varargout = {dims(1), dims(2)};
            end
        end

        function [x, D] = diagMult(v, x, D)
            if any(v)
                if x.rowMajor
                    if x.useMex
                        x.diagonal = mexDiagMult(x.diagonal, v, true);
                    elseif x.allowImplicitExpansion
                        x.diagonal = x.diagonal.*v';
                    else
                        x.diagonal = bsxfun(@times, x.diagonal, v');
                    end
                else
                    if x.allowImplicitExpansion
                        x.diagonal = x.diagonal.*v;
                    else
                        x.diagonal = bsxfun(@times, x.diagonal, v);
                    end
                end
            else
                x = x.toZero();
            end
        end
        
        function [x, D1, D2] = diagProductMult(v1, v2, x, y, D1, D2)
            if isnumeric(x) || isnumeric(y)
                [x, D2] = diagMult(v2, x, D2);
                [y, D1] = diagMult(v1, y, D1);
                x = x + y;
                % Early return
                return
            else
                if isempty(x.diagonal)
                    [x, D1] = diagMult(v1, y, D1);
                elseif isempty(y.diagonal)
                    [x, D2] = diagMult(v2, x, D2);
                else
                    if x.rowMajor
                        if x.useMex
                            x.diagonal = mexDiagProductMult(x.diagonal, v2, y.diagonal, v1, true);
                        elseif x.allowImplicitExpansion
                            x.diagonal = x.diagonal.*v2' + y.diagonal.*v1';
                        else
                            x.diagonal = bsxfun(@times, x.diagonal, v2') + bsxfun(@times, y.diagonal, v1');
                        end
                    else
                        if x.allowImplicitExpansion
                            % Both are diagonal, new Matlab
                            x.diagonal = x.diagonal.*v2 + y.diagonal.*v1;
                        else
                            % Both are diagonal, old Matlab
                            x.diagonal = bsxfun(@times, x.diagonal, v2) + bsxfun(@times, y.diagonal, v1);
                        end
                    end
                end
            end
            sx = x.subset;
            sy = y.subset;
            if isempty(sy)
                x.subset = [];
            elseif ~isempty(sx)
                % Remove zero entries not present in either.
                x.subset = max(sx, sy);
            end
        end
        
        function x = sum(D, n)
            if nargin == 1
                n = 1;
            end
            if size(D.diagonal, 1+D.rowMajor) == 1 && n == 1
                x = D;
            else
                x = sum(D.sparse(), n);
            end
        end
        
        function n = nnz(D)
            n = nnz(D.diagonal);
        end
        
        function u = horzcat(varargin)
            u = DiagonalJacobian.cat(2, varargin{:});
        end
        function u = vertcat(varargin)
            u = DiagonalJacobian.cat(1, varargin{:});
        end
    end
    
    methods(Static)
        function out = cat(dim, varargin)
            if dim == 2 || ...
                    any(cellfun(@isnumeric, varargin)) ||...
                    any(cellfun(@(x) isa(x, 'FixedWidthJacobian'), varargin))
                nz = zeros(1, numel(varargin));
                for i = 1:numel(varargin)
                    nz(i) = nnz(varargin{i});
                end
                cz = cumsum([0, nz]);
                [I, J, V] = deal(zeros(sum(nz), 1));
                [N, M] = deal(0);
                for arg = 1:numel(varargin)
                    ix = (cz(arg)+1):cz(arg+1);
                    if dim == 1
                        [I(ix), J(ix), V(ix), n, M] = getSparseArguments(varargin{arg}, N);
                        N = N + n;
                    else
                        [I(ix), J(ix), V(ix), N, m] = getSparseArguments(varargin{arg}, 0, M);
                        M = M + m;
                    end
                end
                out = sparse(I, J, V, N, M);
            else
                m = numel(varargin);
                l = varargin{1}.dim(2);
                assert(~varargin{1}.rowMajor);
                subsets = cellfun(@(x) x.getSubset(), varargin, 'UniformOutput', false);
                diags = cellfun(@(x) x.diagonal, varargin, 'UniformOutput', false);
                isZero = cellfun(@(x) x.isZero, varargin);
                for i = 1:m
                    if isZero(i)
                        nz = size(diags{i}, 1);
                        diags{i} = zeros(nz, l);
                    end
                end

                d = vertcat(diags{:});
                map = vertcat(subsets{:});
                
                out = varargin{1};
                out.diagonal = d;
                out.subset = map;
            end
        end
    end
    
    methods(Static, Access=private)
        function eq = compareIndices(u, ind_into_u, b_subset)
            if isempty(u.subset)
                % The subset is equal to the indexing
                u_subset = ind_into_u;
            else
                % We get the index of the subset itself
                u_subset = u.getSubset();
                u_subset = u_subset(ind_into_u);
            end
            if islogical(u_subset) && islogical(b_subset)
                eq = all(u_subset(:) == b_subset(:));
                return
            elseif islogical(b_subset)
                b_subset = find(b_subset);
            elseif islogical(u_subset)
                u_subset = find(u_subset);
            end
            eq = all(u_subset(:) == b_subset(:) | u_subset(:) == 0 | b_subset(:) == 0);
        end

        function isZ = isAllZeros(v)
            isZ = isnumeric(v) && ~any(any(v));
        end

        function xsubset = treatSubset(xsubset, ysubset)
            if ~isempty(xsubset) && ~isempty(ysubset)
                isWildcard = xsubset == 0;
                xsubset(isWildcard) = ysubset(isWildcard);
            end
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
