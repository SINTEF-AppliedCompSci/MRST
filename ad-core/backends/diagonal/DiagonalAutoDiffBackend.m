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
        modifyOperators = true; % Update the operators and use custom versions that return `FixedWidthJacobian` instances
        useMex = mrstSettings('get', 'useMEX');
        rowMajor = false;
        deferredAssembly = false;
    end
    
    methods
        function backend = DiagonalAutoDiffBackend(varargin)
            backend = backend@AutoDiffBackend();
            backend = merge_options(backend, varargin{:});
        end
        
        function model = updateDiscreteOperators(backend, model)
            ops = model.operators;
            if isempty(ops) || ~backend.modifyOperators
                % No need to modify empty operators
                return
            end
            if ~isfield(ops, 'N')
                % Need neighborship for this to work
                return
            end
            flds = {'useMex', 'deferredAssembly', 'rowMajor'};
            if isfield(ops, 'operator_parameters')
                operator_parameters = ops.operator_parameters;
                ok = true;
                for i = 1:numel(flds)
                    f = flds{i};
                    if backend.(f) ~= operator_parameters.(f)
                        ok = false;
                        break
                    end
                end
                if ok
                    % Operators match settings, we return
                    return
                end
            end
            N = ops.N;
            nc = model.G.cells.num;
            nf = size(N, 1);
            n1 = N(:, 1);
            n2 = N(:, 2);
            % Discrete divergence
            if isfield(ops, 'Div') || isfield(ops, 'AccDiv')
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

                ops.Div = @(v) discreteDivergence([], v, div_options);
                if numel(N)
                    adiv = @(a, v) discreteDivergence(a, v, div_options);
                else
                    adiv = @(a, v) a;
                end
                ops.AccDiv = adiv;
            end
            % Cell -> Face operators: Grad, upstream and face average
            if isfield(ops, 'Grad')
                ops.Grad = @(v) twoPointGradient(N, v, gradMat, backend.useMex);
            end
            if isfield(ops, 'faceUpstr')
                ops.faceUpstr = @(flag, v) singlePointUpwind(flag, N, v, backend.useMex);
            end
            if isfield(ops, 'faceAvg')
                ops.faceAvg = @(v) faceAverage(N, v, backend.useMex);
            end
            % Store backend operators settings to be somewhat robust if
            % someone modifies them and calls this function again.
            up = struct();
            for i = 1:numel(flds)
                f = flds{i};
                up.(f) = backend.(f);
            end
            ops.operator_parameters = up;
            model.operators = ops;
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
            if ~backend.modifyOperators
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
        
        function [A_cc, b_c, A_nn, b_n, A_cn, A_nc] = getBlockSystemCSR(backend, problem, model, bz)
            if nargin < 3 || bz == 0
                bz = problem.countOfType('cell');
            end
            eqs = problem.equations;
            ne = numel(eqs);
            nj = numel(eqs{1}.jac);
            c_sub = 1:bz; % Assuming cell equations come first
            cell_eq = eqs(c_sub);
            J = applyFunction(@(x) x.jac{1}, cell_eq);
            
            % First, get the matrix
            opt = J{1}.divergenceOptions;
            prelim = opt.mex;
            acc = applyFunction(@(x) x.accumulation.diagonal, J);
            flux = applyFunction(@(x) x.flux.diagonal, J);
            if ~backend.rowMajor
                % MEX routine for assembly assumes rowmajor. Perform
                % potentially expensive transpose of all block matrices.
                acc = cellfun(@(x) x', acc, 'UniformOutput', false);
                flux = cellfun(@(x) x', flux, 'UniformOutput', false);
            end
            bad = cellfun(@isempty, flux);
            if any(bad)
                tmp = flux{bad};
                [flux{bad}] = deal(zeros(2*bz, size(tmp, 2)));
            end
            
            [colNo, rowPtr, values, n, m] = mexDiscreteDivergenceBlockJac(acc, flux, opt.N,...
            prelim.facePos, prelim.faces, prelim.cells, prelim.cellIndex, bz, backend.rowMajor);
            % Get the right hand side
            b_c = getRHS(cell_eq);
            assert(n == m)
            A_cc = struct('col_no', colNo, 'row_ptr', rowPtr, 'val', values, 'n', n, 'block_size', bz, 'mapping', prelim);
            if nargout > 2
                if ne == bz
                    % No extra equations, but they were requested. Just
                    % output empty arrays.
                    [A_nn, b_n, A_cn, A_nc] = deal([]);
                else
                    n_sub = (bz+1):ne;
                    % Assume system:
                    % A = [A_cc, A_cn
                    %      A_nc, A_nn]
                    % where c are cell equations/variables and n are non-cell.
                    n_eq = eqs(n_sub);

                    if nj > 1
                        A_nn = getMatrix(n_eq, 2:nj);
                        if nargout > 3
                            b_n = getRHS(n_eq);
                            if nargout > 4
                                A_cn = getMatrix(cell_eq, 2:nj);
                                if nargout > 5
                                    A_nc = getMatrix(n_eq, 1);
                                end
                            end
                        end
                    else
                        % "Slim" matrix, just two blocks
                        [A_nn, b_n, A_cn] = deal([]);
                        if nargout > 5
                            A_nc = getMatrix(n_eq, 1);
                        end
                    end
                end
            end
        end
        
        function [A, b, schur_diag, insert, fill] = applySchurComplementBlockSystemCSR(backend, A, b, A_nn, A_nc, A_cn, b_n, schur_type, w)
            if nargin < 9
                w = 1;
            end
            if nargin < 8
                doRow = true;
                doSum = false;
            else
                switch schur_type
                    case 'diagonal'
                        doRow = true;
                        doSum = false;
                    case 'rowsum'
                        doRow = true;
                        doSum = true;
                    case 'colsum'
                        doRow = false;
                        doSum = true;
                    otherwise
                        error('% is not supported.', schur_type);
                end
            end
            bz = A.block_size;
            n = size(A_nc, 2)/bz;
            if isnumeric(A_nn)
                n_to_c = A_nn\A_nc;
                x_n = A_nn\b_n;
            else
                n_to_c = A_nn(A_nc);
                x_n = A_nn(b_n);
            end
            b = b - A_cn*x_n;
            fill = A_cn*n_to_c;
            [row, col, vals] = find(fill);
            cell_row = mod(row-1, n) + 1;
            cell_col = mod(col-1, n) + 1;
            keep = cell_row == cell_col;
            ix = unique(cell_row(keep));
            if doSum
                % vals(~keep) = -vals(~keep);
            else
                row = row(keep);
                col = col(keep);
                vals = vals(keep);
            end
            schur_diag = zeros(bz*bz, numel(ix));

            % Block
            block_row = floor((row-1)/n) + 1;
            block_col = floor((col-1)/n) + 1;
            % Create block values to insert
            for eqNo = 1:bz
                for derNo = 1:bz
                    act = block_row == eqNo & block_col == derNo;
                    if doRow
                        ii = row(act) - (eqNo-1)*n;
                    else
                        ii = col(act) - (derNo-1)*n;
                    end
                    rowsum = accumarray(ii, vals(act), [n, 1]);
                    sub = (eqNo-1)*bz + derNo;
                    schur_diag(sub, :) = schur_diag(sub, :) + rowsum(ix)';
                end
            end
            prelim = A.mapping;
            d_offset = prelim.cellIndex(ix);
            insert = prelim.facePos(ix) + d_offset + ix;
            if w ~= 0
                A.val(:, insert) = A.val(:, insert) - w.*schur_diag;
            end
        end
    end
end



function A = getMatrix(eqs, j_ix)   
    A = cell(1, numel(eqs));
    for i = 1:numel(eqs)
        A{i} = horzcat(eqs{i}.jac{j_ix});
    end
    A = vertcat(A{:});
end

function b = getRHS(eqs)
    values = cellfun(@(x) -value(x), eqs, 'UniformOutput', false);
    b = vertcat(values{:});
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
