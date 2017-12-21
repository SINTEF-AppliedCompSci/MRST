classdef DiagonalSubset < DiagonalJacobian
    % Structured subset of a diagonal jacobian 
    properties
        map % Map to the underlying DiagonalJacobian representation. Two DiagonalSubsets of the same map can be multiplied together, etc.
    end
    
    methods
        function D = DiagonalSubset(d, dims, map, subset)
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
            end
        end

        function [x, D] = diagMult(v, x, D)
            if any(x.diagonal(:))
                v = repmat(v, size(x.map, 2), 1);
                x.diagonal = bsxfun(@times, x.diagonal, v);
            else
                x = x.toZero();
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
            imax = n/nmap;
            jmax = prod(D.dim);
            if isempty(D.subset)
                nval = size(D.diagonal, 1);
            else
                nval = numel(D.subset)/nmap;
            end
            if size(D.map, 1) == 1
                D.map = repmat(D.map, nval, 1);
            end
            I = repmat((1:n/nmap)', nmap, m);
            if isempty(D.subset)
                jmap = reshape(D.map, [], 1);
            else
                jmap = reshape(D.map(D.subset, :), [], 1);
            end
            
            if m == 1
                J = reshape(jmap, [], 1);
            else
                J = bsxfun(@plus, repmat(jmap, 1, m), (0:m-1)*D.dim(1));
            end
            V = D.diagonal;
        end
        
        function u = repmat(u, varargin)
            u = repmat(u.sparse(), varargin{:});
        end
    end
end