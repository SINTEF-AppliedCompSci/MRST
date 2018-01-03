classdef (InferiorClasses = {?DiagonalJacobian,?DiagonalSubset}) DivergenceTerm
    % Very experimental divergence term
    properties
        I
        J
        V
        n
        m
        other
    end
    
    methods
        function D = DivergenceTerm(I, J, V, N, M)
            D.I = I; D.J = J; D.V = V; D.n = N; D.m = M;
        end
        
        function [I, J, V, n, m] = getSparseBlocks(D, varargin)
            if isempty(D.other)
                I = D.I;
                J = D.J;
                V = D.V;
            else
                [Io, Jo, Vo] = getSparseBlocks(D.other);
                I = [D.I(:); Io(:)];
                J = [D.J(:); Jo(:)];
                V = [D.V(:); Vo(:)];
            end
            n = D.n;
            m = D.m;
            if nargin > 1
                I = I + varargin{1};
                if nargin > 2
                    J = J + varargin{2};
                end
            end
        end

        
        function [I, J, V, n, m] = getSparseArguments(D, varargin)
            [I, J, V, n, m] = getSparseBlocks(D);
            
            act = V ~= 0;
            I = I(act);
            J = J(act);
            V = V(act);
            if nargin > 1
                I = I + varargin{1};
                if nargin > 2
                    J = J + varargin{2};
                end
            end
        end
        
        function s = sparse(D)
            [i, j, v, nn, mm] = D.getSparseArguments();
            s = sparse(i, j, v, nn, mm);
        end
        
        function D = plus(D, v)
            if isa(D, 'DivergenceTerm')
                D.other = v;
            else
                v.other = D;
                D = v;
            end
        end
        
        function v = nnz(D)
            v = nnz(D.other) + nnz(D.V);
        end
        
        function [x, D] = diagMult(v, x, D)
            x.V = x.V.*v(x.I);
            if ~isempty(x.other)
                [x.other, D] = diagMult(v, x.other, D);
            end
        end
        
        function D = uminus(D)
            D.V = -D.V;
            D.other = -D.other;
        end
        
        function varargout = matrixDims(D, n)
            dims = [D.n, D.m];
            if nargout == 1
                varargout{1} = dims;
                if nargin > 1
                    varargout{1} = varargout{1}(n);
                end
            else
                varargout = {dims(1), dims(2)};
            end
        end
    end
end