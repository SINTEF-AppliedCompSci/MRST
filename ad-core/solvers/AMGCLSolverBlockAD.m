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
            bz = solver.amgcl_setup.block_size;
            if bz == 0
                % Defaulted - we estimate
                bz = problem.countOfType('cell');
                solver.amgcl_setup.block_size = bz;
            end
            backend = model.AutoDiffBackend;
            [A, b, A_nn, b_n, A_cn, A_nc] = backend.getBlockSystemCSR(problem, model, bz);
            t_asm = toc(timer);
            needsReduction = ~isempty(A_nn);
            if needsReduction
                % We have non-cell equations. We try to do something about
                % it.
                s = solver.reductionStrategy;
                solveBoth = strcmpi(s, 'prepost');
                solvePre  = solveBoth || strcmp(s, 'pre');
                solvePost = solveBoth || strcmp(s, 'post');
                solveSchur = strcmpi(s, 'schur');
                
                assert(solvePre || solvePost || solveSchur);
                % Factor remaining equations
                [L, U] = lu(A_nn);
                A_nn_inv = @(x) U\(L\x);
                if solvePre
                    x_n = A_nn_inv(b_n);
                    b = b - A_cn*x_n;
                elseif solveSchur
                    [A, b] = backend.applySchurComplementBlockSystemCSR(A, b, A_nn_inv, A_nc, A_cn, b_n, ...
                                                                   solver.schurApproxType, solver.schurWeight);
                end
            end
            t_prep = toc(timer);
            % Pass of specialized matrix to MEX file
            [x_c, report] = solver.callAMGCL_MEX(A.col_no, A.row_ptr, A.val, b);
            t_solve = toc(timer) - t_prep;
            if needsReduction
                if solvePost || solveSchur
                    x_n = A_nn_inv(b_n - A_nc*x_c);
                end
                result = [x_c; x_n];
            else
                result = x_c;
            end
            t_post = toc(timer) - t_solve - t_prep;
            % Output helpful info
            report.PreparationTime    = t_prep;
            report.LinearSolutionTime = t_solve;
            report.PostProcessTime    = t_post;
            report.BlockAssembly      = t_asm;
            dx = solver.storeIncrements(problem, result);
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
