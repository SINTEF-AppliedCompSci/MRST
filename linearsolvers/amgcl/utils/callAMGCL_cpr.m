function [x, err, nIter] = callAMGCL_cpr(A, b, block_size, varargin)
%Invoke AMGCL Linear Solver Software in CPR Mode
%
% DESCRIPTION:
%   Solves the system of simultaneous linear equations
%
%      Ax = b
%
%   in which the coefficient matrix 'A' is sparse.  Single right-hand side
%   vector.
%
% SYNOPSIS:
%    x              = callAMGCL_cpr(A, b, block_size)
%    x              = callAMGCL_cpr(A, b, block_size, 'pn1', pv1, ...)
%   [x, err]        = callAMGCL_cpr(...)
%   [x, err, nIter] = callAMGCL_cpr(...)
%
% PARAMETERS:
%   A          - Coefficient matrix.
%
%   b          - System right-hand side vector.
%
%   block_size - Size of sub-blocks.  Positive integer.
%
% KEYWORD ARGUMENTS:
%   isTransposed   - Whether or not the coefficient matrix is transposed on
%                    input.  This simplifies converting MATLAB's CSC matrix
%                    format to AMGCL's CSR matrix input format.  Default
%                    value: `isTransposed = false`.
%
%   cellMajorOrder - Whether or not the coefficient matrix already is in
%                    cell-major ordering (i.e., whether or not the degrees
%                    of freedom are grouped by cell).  If set, this makes
%                    the call into AMGCL proper less expensive than
%                    otherwise.
%
%   tolerance      - Linear solver tolerance.  Corresponds to parameter
%                    'solver.tol' (relative residual reduction) of AMGCL.
%                    Default value: `tolerance = 1.0e-6`.
%
%   maxIterations  - Maximum number of linear iterations.  Overrides
%                    'solver.maxiter' if positive.  Default value:
%                    `maxIterations = 0` (use AMGCL solver's default,
%                    typically 100).
%
%  Additional keyword arguments passed on to function `getAMGCLMexStruct`.
%
% RETURNS:
%   x     - Solution vector.
%
%   err   - Norm of residual at end of solution process.
%
%   nIter - Number of linear iterations.
%
% SEE ALSO:
%   `callAMGCL`, `amgcl_matlab`, `getAMGCLMexStruct`.

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

    assert(block_size > 0);

    opt = struct('isTransposed',   false, ...
                 'cellMajorOrder', false, ...
                 'tolerance',     1e-6, ...
                 'reuseMode',     1, ...
                 'maxIterations',  0);
    [opt, cl_opts] = merge_options(opt, varargin{:});
        
    amg_opt = getAMGCLMexStruct(cl_opts{:});
    amg_opt.block_size = block_size;
    
    assert(mod(size(A, 1), block_size) == 0);
    if opt.cellMajorOrder
        if ~opt.isTransposed
            A = A';
        end
    else
        n = size(A, 1);
        ordering = getCellMajorReordering(n/block_size, block_size, 'ndof', n);
        A = A(ordering, ordering)';
        b = b(ordering);
    end

    t = tic();
    szb = size(b);
    if szb(2) > 1
        % Multiple right-hand-sides
        x = zeros(szb);
        if opt.reuseMode == 1
            resetAMGCL();
        end
        [x(:, 1), err, nIter] = ...
           amgcl_matlab(A, b(:, 1), amg_opt, opt.tolerance, opt.maxIterations, 2, 2);
        for i = 2:szb(2)
            [x(:, i), ~, it] = ...
               amgcl_matlab(A, b(:, i), amg_opt, opt.tolerance, opt.maxIterations, 1, 2);
           nIter = nIter + it;
        end
        nIter = nIter/szb(2);
        if opt.reuseMode == 1
            % Clean up
            resetAMGCL();
        end
    else
        % A single right-hand-side
        [x, err, nIter] = ...
           amgcl_matlab(A, b, amg_opt, opt.tolerance, opt.maxIterations, 2, opt.reuseMode);
    end
    t_solve = toc(t);
    if ~opt.cellMajorOrder
        x(ordering) = x;
    end
    
    if err > opt.tolerance
        warning('Solver did not converge to specified tolerance of %e in %d iterations. Reported residual estimate was %e after %2.2f seconds', opt.tolerance, nIter, err, t_solve);
    elseif mrstVerbose()
        fprintf('AMGCL solver converged to %e in %2d iterations after %2.2f seconds.\n', err, nIter, t_solve);
    end
end
