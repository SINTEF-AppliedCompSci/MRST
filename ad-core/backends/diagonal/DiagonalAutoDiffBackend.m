classdef DiagonalAutoDiffBackend < AutoDiffBackend
    properties
        modifyOperators = true;
    end
    
    methods
        function backend = DiagonalAutoDiffBackend()
            backend = backend@AutoDiffBackend();
        end
        
        function model = updateDiscreteOperators(backend, model)
            if isa(model, 'ReservoirModel') && backend.modifyOperators
                N = model.operators.N;
                nc = model.G.cells.num;
                nf = size(N, 1);
                n1 = N(:, 1);
                n2 = N(:, 2);


                C  = sparse(N, [(1:nf)'; (1:nf)'], ones(nf,1)*[1 -1], nc, nf);

                gradMat = sparse((1:2*nf), [n1; n2], rldecode([-1; 1], [nf; nf]), 2*nf, nc);

                [~, sortedN] = sort(repmat(reshape(N, [], 1), 2, 1));
                I_base = [N(:, 1); N(:, 1); N(:, 2); N(:, 2)];
                I_base = I_base(sortedN);

                sortIx = struct('C', C, 'J_sorted_index', sortedN, 'I_base', I_base);
                model.operators.Grad = @(v) twoPointGradient(N, v, gradMat);
                model.operators.Div = @(v) discreteDivergence(N, v, nc, nf, sortIx, C);
                model.operators.faceUpstr = @(flag, v) singlePointUpwind(flag, N, v);
                model.operators.faceAvg = @(v) faceAverage(N, v);
            end
        end
        
        function out = getBackendDescription(backend)
            out = 'Diagonal AD';
            if backend.modifyOperators
                out = [out, ' [custom operators]'];
            else
                out = [out, ' [matrix product operators]'];
            end
        end
        function varargout = initVariablesAD(backend, varargin)
           n         = nargout;
           varargout = cell([1, n]);
           
           [varargout{:}] = initVariablesAD_diagonal(varargin{:});
        end
        
        function v = convertToAD(backend, v, sample)
           v = double2NewAD(v, sample);
        end
    end
end