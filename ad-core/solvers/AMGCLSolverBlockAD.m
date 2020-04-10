classdef AMGCLSolverBlockAD < AMGCLSolverAD
    % Linear solver that calls external compiled multigrid solver
    %
    % SYNOPSIS:
    %   solver = AMGCLSolverAD()
    %
    % DESCRIPTION:
    %    AD-interface for the AMGCL interface.
    %
    % NOTE:
    %    This solver requires AMGCL to be installed and working.
    %
    % SEE ALSO:
    %   `BackslashSolverAD`
    properties
        reductionStrategy = 'schur';
        schurApproxType = 'diagonal';
        schurWeight = 1;
    end
    
    properties (Access = protected)
        amgcl_id = 1;
    end
    
    methods
        function solver = AMGCLSolverBlockAD(varargin)
            solver = solver@AMGCLSolverAD(varargin{:});
        end
        
        function [result, report] = solveLinearSystem(solver, A, b)
            assert(false);
        end
        
        function [dx, result, report] = solveLinearProblem(solver, problem, model)
            timer = tic();
            backend = model.AutoDiffBackend;
            eqs = problem.equations;
            bz = solver.amgcl_setup.block_size;
            if bz == 0
                % Defaulted - we estimate
                bz = problem.countOfType('cell');
                solver.amgcl_setup.block_size = bz;
            end
            ne = numel(eqs);
            nj = numel(eqs{1}.jac);
            needsReduction = ne > bz;
            c_sub = 1:bz; % Assuming cell equations come first
            cell_eq = eqs(c_sub);
            J = cellfun(@(x) x.jac{1}, cell_eq, 'UniformOutput', false);
            
            % First, get the matrix
            opt = J{1}.divergenceOptions;
            prelim = opt.mex;
            
            acc = cellfun(@(x) x.accumulation.diagonal, J, 'UniformOutput', false);
            flux = cellfun(@(x) x.flux.diagonal, J, 'UniformOutput', false);
            bad = cellfun(@isempty, flux);
            if any(bad)
                [flux{bad}] = deal(zeros(2*bz, size(model.operators.N, 1)));
            end
            
            [colNo, rowPtr, values, n, m] = mexDiscreteDivergenceBlockJac(acc, flux, opt.N,...
            prelim.facePos, prelim.faces, prelim.cells, prelim.cellIndex, bz, backend.rowMajor);
            % Get the right hand side
            b = getRHS(cell_eq);
            assert(n == m)
            if needsReduction
                % We have non-cell equations. We try to do something about
                % it.
                s = solver.reductionStrategy;
                solveBoth = strcmpi(s, 'prepost');
                solvePre  = solveBoth || strcmp(s, 'pre');
                solvePost = solveBoth || strcmp(s, 'post');
                solveSchur = strcmpi(s, 'schur');
                
                assert(solvePre || solvePost || solveSchur);
                n_sub = (bz+1):ne;
                % Assume system:
                % A = [A_cc, A_cn
                %      A_nc, A_nn]
                % where c are cell equations/variables and n are non-cell.
                n_eq = eqs(n_sub);
                A_cn = getMatrix(cell_eq, 2:nj);
                A_nn = getMatrix(n_eq, 2:nj);
                if solvePost || solveSchur
                    A_nc = getMatrix(n_eq, 1);
                end
                b_n = getRHS(n_eq);
                % Factor remaining equations
                [L, U] = lu(A_nn);
                A_nn_inv = @(x) U\(L\x);
                if solvePre || solveSchur
                    % A_nn_inv = @(x) A_nn\x;
                    x_n = A_nn_inv(b_n);
                    b = b - A_cn*x_n;
                    if solveSchur
                        tic()
                        [schur_diag, ix] = solver.getSchurBlockComplement(A_nn_inv, A_nc, A_cn, bz);
                        d_offset = prelim.cellIndex(ix);
                        insert = prelim.facePos(ix) + d_offset + ix;
                        values(:, insert) = values(:, insert) - schur_diag;
                        toc()
                    end
                end
            end
            t_prep = toc(timer);
            % Pass of specialized matrix to MEX file
            [x_c, report] = solver.callAMGCL_MEX(colNo, rowPtr, values, b);
            t_solve = toc(timer) - t_prep;
            if needsReduction
                if solvePost || solveSchur
                    x_n = A_nn_inv(b_n - A_nc*x_c);
                end
                result = [x_c; x_n];
            else
                result = x_c;
            end
            dx = solver.storeIncrements(problem, result);
            t_post = toc(timer) - t_solve;
            % Output helpful info
            report.PreparationTime = t_prep;
            report.LinearSolutionTime = t_solve;
            report.PostProcessTime = t_post;
        end
        
        
        function [result, report] = callAMGCL_MEX(solver, I, J, V, b)
            tol = solver.tolerance;
            timer = tic();
            [result, res, its] = amgcl_matlab_block(I, J, V, b, solver.amgcl_setup, ...
                                                    tol, solver.maxIterations, solver.amgcl_id);
            t_solve = toc(timer);
            if res > solver.tolerance
                warning(['Solver did not converge to specified tolerance of %1.3e in %d iterations. ', ...
                    'Reported residual estimate was %1.3e after %2.2f seconds'], tol, its, res, t_solve);
            elseif solver.verbose
                fprintf('AMGCL solver converged to %1.3e in %2d iterations after %2.2f seconds.\n', res, its, t_solve);
            end
            report = solver.getSolveReport(...
                'Converged',  res <= solver.tolerance, ...
                'Residual',   res,...
                'Iterations', its);
        end
        
        function [schur_diag, ix] = getSchurBlockComplement(solver, A_nn_inv, A_nc, A_cn, bz)
            st = solver.schurApproxType;
            switch st
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
                    error('% is not supported.', st);
            end
            n = size(A_nc, 2)/bz;
            n_to_c = A_nn_inv(A_nc);
            fill = A_cn*n_to_c;
            [row, col, vals] = find(fill);
            cell_row = mod(row-1, n) + 1;
            cell_col = mod(col-1, n) + 1;
            keep = cell_row == cell_col;
            ix = unique(cell_row(keep));
            if ~doSum
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
                    sub = (derNo-1)*bz + eqNo;
                    schur_diag(sub, :) = schur_diag(sub, :) + rowsum(ix)';
                end
            end
            w = solver.schurWeight; % Pessimist factor for approximate schur
            if w ~= 0
                schur_diag = schur_diag.*w;
            end
            % l = 6; k = ix(l); full(fill(k + [0, 1, 2]*n, k + [0, 1, 2]*n)), schur_diag(:, l)
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
