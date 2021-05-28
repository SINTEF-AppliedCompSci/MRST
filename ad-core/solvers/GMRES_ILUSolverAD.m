classdef GMRES_ILUSolverAD < LinearSolverAD
    %Preconditioned GMRES solver.
    %
    % SYNOPSIS:
    %   solver = GMRES_ILUSolverAD()
    %
    % DESCRIPTION:
    %   Solve a linearized problem using a GMRES solver with ILU
    %   preconditioner.
    %
    % REQUIRED PARAMETERS:
    %   None
    %
    % OPTIONAL PARAMETERS:
    %   'property' - Set property to value.
    %
    %
    % SEE ALSO:
    %   `BackslashSolverAD`, `CPRSolverAD`, `LinearizedProblem`,
    %   `LinearSolverAD`

    properties
        % Type: {'nofill'}, 'crout', 'ilutp'
        ilutype
        % Drop tolerance for elements (default: 0)
        dropTolerance
        % Modified ilu type: 'row', 'col', {'off'}
        modifiedIncompleteILU
        % Boolean indicating replacement of zero diagonal entries
        udiagReplacement
        % Pivoting threshold for ilupt. Default 1.
        pivotThreshold
        % Reorder equations when diagonal entries are zero for ilu0
        reorderEquations
    end
    methods
        function solver = GMRES_ILUSolverAD(varargin)
            solver = solver@LinearSolverAD();
            
            solver.ilutype               = 'nofill';
            solver.dropTolerance         = 0;
            solver.modifiedIncompleteILU = 'off';
            solver.udiagReplacement      = true;
            solver.pivotThreshold        = 1;
            solver.reorderEquations      = true;
            
            solver = merge_options(solver, varargin{:});
        end
        
        function [result, report] = solveLinearSystem(solver, A, b, varargin)
            timer = tic();
            nel = size(A, 1);
            if solver.reorderEquations
                [A, b] = reorderForILU(A, b);
            end
            [L, U] = ilu(A, solver.getOptsILU());
            [result, flag, res, its] = gmres(A, b, [], ...
                solver.tolerance,...
                min(solver.maxIterations, nel), ...
                L, U, ...
                varargin{:});
            t_solve = toc(timer);
            report = solver.getSolveReport('Residual',   res,...
                                           'Iterations', its, ...
                                           'Converged', flag == 0, ...
                                           'LinearSolutionTime', t_solve);
            report.GMRESFlag = flag;
            if flag > 0
                warning('Solver did not converge to specified tolerance of %g. Reported residual estimate was %g', solver.tolerance, res);
            end
        end
        
        function opts = getOptsILU(solver)
            opts = struct('type',    solver.ilutype, ...
                'droptol', solver.dropTolerance, ...
                'milu',    solver.modifiedIncompleteILU, ...
                'udiag',   solver.udiagReplacement, ...
                'thresh',  solver.pivotThreshold);
        end
        
        function [d, sn] = getDescription(solver)
            sn = 'Matlab-GMRES-ILU(0)';
            sn = [sn, solver.id];
            d = 'GMRES with incomplete LU-factorization (ILU(0)) as preconditioner';
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
