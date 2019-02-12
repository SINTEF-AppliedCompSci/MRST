classdef DiagonalAutoDiffBackend < AutoDiffBackend
    % Automatic differentiation backend class (diagonal representation)
    %
    % SYNOPSIS:
    %   backend = DiagonalAutoDiffBackend()
    %
    % DESCRIPTION:
    %   This backend uses a diagonal representation (with an optional set
    %   of operators that produce intermediate diagonal representations).
    %   The primary use of this class is for problems with a large number
    %   of independent primary variables, with many cell-wise operations
    %   (e.g. compositional or similar problems).
    %
    % RETURNS:
    %   Backend - Initialized class instance
    %
    % SEE ALSO:
    %   `AutoDiffBackend`, `SparseAutoDiffBackend`

    properties
        modifyOperators = true; % Update the operators and use custom versions that return `DiagonalSubset` instances
        useMex = false;
    end
    
    methods
        function backend = DiagonalAutoDiffBackend(varargin)
            backend = backend@AutoDiffBackend();
            backend = merge_options(backend, varargin{:});
        end
        
        function model = updateDiscreteOperators(backend, model)
            if isa(model, 'ReservoirModel') && backend.modifyOperators ...
                    && ~isfield(model.operators, 'diag_updated')
                N = model.operators.N;
                nc = model.G.cells.num;
                nf = size(N, 1);
                n1 = N(:, 1);
                n2 = N(:, 2);


                C  = sparse(N, [(1:nf)'; (1:nf)'], ones(nf,1)*[1 -1], nc, nf);
                gradMat = sparse((1:2*nf), [n1; n2], rldecode([-1; 1], [nf; nf]), 2*nf, nc);
                [~, sortedN] = sort(repmat(reshape(N, [], 1), 2, 1));
                I_base = [N(:, 1); N(:, 1); N(:, 2); N(:, 2)];
                if backend.useMex
                    prelim = getMexDiscreteDivergenceJacPrecomputes(model);
                else
                    prelim = [];
                end

                sortIx = struct('C', C, 'J_sorted_index', sortedN, 'I_base', I_base);
                model.operators.Grad = @(v) twoPointGradient(N, v, gradMat, backend.useMex);
                model.operators.Div = @(v) discreteDivergence(N, v, nc, nf, sortIx, C, prelim, backend.useMex);
                model.operators.faceUpstr = @(flag, v) singlePointUpwind(flag, N, v, backend.useMex);
                model.operators.faceAvg = @(v) faceAverage(N, v, backend.useMex);
                model.operators.diag_updated = true;
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
           if nargout ~= nargin
               opts = varargin{end};
               varargin = varargin(1:end-1);
           else
               opts = struct('useMex', backend.useMex, 'types', []);
           end
           [varargout{:}] = initVariablesAD_diagonal(varargin{:}, opts);
        end
        
        function v = convertToAD(backend, v, sample)
           v = double2NewAD(v, sample);
        end
    end
end
