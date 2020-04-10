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
        rowMajor = false;
        deferredAssembly = false;
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
                % Discrete divergence
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
                div_options = struct('sortIx',             sortIx,...
                                     'mex',                prelim, ...
                                     'useMex',             backend.useMex, ...
                                     'fn',                 [], ...
                                     'useConservationJac', backend.deferredAssembly, ...
                                     'N', N, 'C', C, 'nf', nf, 'nc', nc);

                model.operators.Div = @(v) discreteDivergence([], v, div_options);
                model.operators.AccDiv = @(a, v) discreteDivergence(a, v, div_options);
                % Cell -> Face operators: Grad, upstream and face average
                model.operators.Grad = @(v) twoPointGradient(N, v, gradMat, backend.useMex);
                model.operators.faceUpstr = @(flag, v) singlePointUpwind(flag, N, v, backend.useMex);
                model.operators.faceAvg = @(v) faceAverage(N, v, backend.useMex);

                model.operators.diag_updated = true;
            end
        end
        
        function out = getBackendDescription(backend)
            if backend.useMex
                out = 'Diagonal-MEX';
            else
                out = 'Diagonal';
            end
            if backend.rowMajor
                out = [out, '-RowMajor'];
            end
            if backend.modifyOperators
                out = [out, ' [custom operators]'];
            else
                out = [out, ' [matrix product operators]'];
            end
        end
        function varargout = initVariablesAD(backend, varargin)
           n         = nargout;
           varargout = cell([1, n]);
           if nargin > nargout + 1
               opts = varargin{end};
               varargin = varargin(1:end-1);
           else
               opts = struct('useMex', backend.useMex, 'types', []);
           end
           if backend.rowMajor
               [varargout{:}] = initVariablesAD_diagonalRowMajor(varargin{:}, opts);
           else
               [varargout{:}] = initVariablesAD_diagonal(varargin{:}, opts);
           end
        end
        
        function v = convertToAD(backend, v, sample)
           v = double2GenericAD(v, sample);
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
