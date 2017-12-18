classdef DiagonalJacobian
    properties
        diagonal
        subset % List of subset values
        dim % Vector: First dimension is the number of variables in block, while the second is the number of columns
    end
    
    methods
        function D = DiagonalJacobian(d, dim, subset)
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
        end
        
        function sub = getSubset(u)
            if u.isZero
                sub = zeros(size(u.diagonal, 1), 1);
            elseif isempty(u.subset)
                sub = reshape(1:size(u.diagonal, 1), [], 1);
            else
                sub = u.subset;
            end
        end

        function u = toZero(u, n)
            if nargin == 1
                n = size(u.diagonal, 1);
            end
            u.diagonal = zeros(n, 0);
            % Unsure if this is the best idea.
            u.subset = zeros(n, 1);
        end

        function u = expandZero(u)
            u.diagonal = zeros(size(u.diagonal, 1), u.dim(2));
        end
        
        function isz = isZero(u)
            isz = size(u.diagonal, 2) == 0;
        end
        
        function [I, J, V, imax, jmax] = getSparseBlocks(D)
            [imax, m] = size(D.diagonal);
            jmax = prod(D.dim);
            if m == 0
                I = []; J = []; V = [];
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
            V = D.diagonal;
        end
        
        function [I, J, V, imax, jmax] = getSparseArguments(D, ioffset, joffset)
            [I, J, V, imax, jmax] = getSparseBlocks(D);
            act = V ~= 0;
            I = I(act);
            J = J(act);
            V = V(act);
            if nargin > 1
                I = I + ioffset;
                if nargin > 2
                    J = J + joffset;
                end
            end
        end
        
        function s = sparse(D)
            [I, J, V, n, m] = D.getSparseArguments();
            s = sparse(I, J, V, n, m);
        end
        
        function s = double(D)
            s = D.sparse();
        end

        function isEqual = subsetsEqualNoZeroCheck(x, y, xsubset, ysubset)
            if nargin == 2
                xsubset = x.subset;
                ysubset = y.subset;
            end
            chx = isempty(xsubset);
            chy = isempty(ysubset);
            if chx && chy
                isEqual = true;
            elseif chx
                tmp = reshape(1:x.dim(1), [], 1);
                isEqual = numel(ysubset) == x.dim(1) && all(y.subset(ysubset ~= 0) == tmp(ysubset ~= 0));
            elseif chy
                tmp = reshape(1:x.dim(1), [], 1);
                isEqual = numel(xsubset) == y.dim(1) && all(xsubset(xsubset ~= 0) == tmp(xsubset ~= 0));
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
            u.diagonal = repmat(u.diagonal, varargin{:});
            u.subset = repmat(u.subset, varargin{:});
        end
        function x = subsetPlus(x, v, subs)
            if isa(x, 'DiagonalJacobian')
                if isa(v, 'DiagonalJacobian')
                    if v.isZero
                        return
                    end
                    if x.isZero
                        x = x.sparse();
                        x(subs, :) = x(subs, :) + v;
                        return
                    end
                    % Get subset, check individual values
                    s = x.getSubset();
                    if subsetsEqualNoZeroCheck(x, v, s(subs), v.subset)
                        x.diagonal(subs, :) = x.diagonal(subs, :) + v.diagonal;
                        if ~isempty(x.subset) && ~isempty(y.subset)
                            isWildcard = x.subset == 0;
                            x.subset(isWildcard) = y.subset(isWildcard);
                        end
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
                        if u.isZero
                            u = u.toZero(numel(s(1).subs{1}));
                        elseif numel(s(1).subs) == 2 && ischar(s(1).subs{2})
                            if any(s(1).subs{1})
                                u.diagonal = u.diagonal(s(1).subs{1}, :);
                            else
                                u = u.toZero(0);
                            end
                            
                            if isempty(u.subset)
                                u.subset = reshape(s(1).subs{1}, [], 1);
                            else
                                subs = u.getSubset();
                                u.subset = subs(s(1).subs{1});
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
                            if u.isZero && v.isZero
                                return
                            end
                            if v.isZero
                                u.diagonal(s.subs{1}, :) = 0;
                                return
                            end
                            if u.isZero
                                u = u.expandZero();
                                u.subset(s.subs{1}) = v.getSubset();
                                u.diagonal(s.subs{1}, :) = v.diagonal;
                                return
                            end
                            if subsetsEqualNoZeroCheck(u, v)
                                u.diagonal(s.subs{1}, :) = v.diagonal;
                            elseif isempty(u.subset) && ~isempty(v.subset) && u.compareIndices(s.subs{1}, v.subset)
                                u.diagonal(s.subs{1}, :) = v.diagonal;
                            else
                                u = u.sparse();
                                v = v.sparse();
                                u(s.subs{1}, :) = v;
                            end
                        else
                            if uD
                                if numel(v) == 1
                                    u.diagonal(s.subs{1}, :) = v;
                                    return
                                end
                                u = u.sparse();
                                % This case can be optimized more.
                            end
                            if vD
                                v = v.sparse();
                            end
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
                    if ~isempty(x.subset)
                        isWildcard = x.subset == 0;
                        x.subset(isWildcard) = y.subset(isWildcard);
                    end
                else
                    x = x.sparse() + y.sparse();
                end
                return
            end
            
            if xD
                if ~DiagonalJacobian.isAllZeros(y)
                    x = x.sparse();
                    x = plus(x, y);
                elseif isscalar(y)
                    x.diagonal = x.diagonal + y;
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
                    if ~isempty(x.subset) && ~isempty(y.subset)
                        isWildcard = x.subset == 0;
                        x.subset(isWildcard) = y.subset(isWildcard);
                    end
                else
                    x = x.sparse().*y.sparse();
                end
                return
            end
            
            if xD
                if isscalar(y)
                    x.diagonal = x.diagonal.*y;
                else
                    x = x.sparse();
                    x = x.*y;
                end
            else
                x = times(y, x);
            end
        end
        
        function x = mtimes(x, y)
            if DiagonalJacobian.isAllZeros(x) || DiagonalJacobian.isAllZeros(y)
                dx = matrixDims(x);
                dy = matrixDims(y);
                assert(dx(2) == dy(1))
                x = sparse([], [], [], dx(1), dy(2));
                return;
            end
            
            xD = isa(x, 'DiagonalJacobian');
            yD = isa(y, 'DiagonalJacobian');
            if xD && yD
                if subsetsEqualNoZeroCheck(x, y)
                    x.diagonal = x.diagonal.*y.diagonal;
                    if ~isempty(x.subset) && ~isempty(y.subset)
                        isWildcard = x.subset == 0;
                        x.subset(isWildcard) = y.subset(isWildcard);
                    end
                else
                    x = x.sparse().*y.sparse();
                end
            elseif xD
                if isscalar(y)
                    x.diagonal = x.diagonal.*y;
                else
                    x = x.sparse();
                    x = x.*y;
                end
            else
                if isscalar(x)
                    y.diagonal = x.*y.diagonal;
                    x = y;
                else
                    y = y.sparse();
                    x = x*y;
                end
            end
        end
        
        function varargout = matrixDims(D, n)
            dims = [size(D.diagonal, 1), prod(D.dim)];
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
            x.diagonal = bsxfun(@times, x.diagonal, v);
        end
        
        function x = sum(D, n)
            if nargin == 1
                n = 1;
            end
            
            if size(D.diagonal, 1) == 1 && n == 1
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
                    any(cellfun(@(x) isa(x, 'DiagonalSubset'), varargin))
                nz = cellfun(@nnz, varargin);
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
                
                subsets = cellfun(@(x) x.getSubset, varargin, 'UniformOutput', false);
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
        function eq = compareIndices(a, b)
            if islogical(a) && islogical(b)
                eq = all(a(:) == b(:));
                return
            elseif islogical(b)
                b = find(b);
            elseif islogical(a)
                a = find(a);
            end
            eq = all(a(:) == b(:));
        end

        function isZ = isAllZeros(v)
            isZ = isnumeric(v) && ~any(v(:));
            %     isZ = nnz(v) == 0;
        end

    end
end

