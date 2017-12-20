classdef (InferiorClasses = {?DiagonalJacobian,?DiagonalSubset}) DivergenceTerm
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
        
        function [I, J, V, n, m] = getSparseArguments(D)
            if isempty(D.other)
                I = D.I;
                J = D.J;
                V = D.V;
            else
                [Io, Jo, Vo] = getSparseArguments(D.other);
                I = [D.I; reshape(Io, [], 1)];
                J = [D.J; reshape(Jo, [], 1)];
                V = [D.V; reshape(Vo, [], 1)];
            end
            n = D.n;
            m = D.m;
        end
        
        function s = sparse(D)
            [i, j, v, nn, mm] = D.getSparseArguments();
            s = sparse(i, j, v, nn, mm);
        end
        
        function D = plus(D, v)
            if isa(D, 'DivergenceTerm')
                if isa(v, 'DiagonalJacobian')
                    D.other = v;
                else
                    D = D.sparse() + v;
                end
            else
                D = plus(v, D);
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