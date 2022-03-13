function [p_ms, report] = solveMultiscaleIteratively(A, q, p_ms, basis, getSmootherFn, tol, iterations, LinSolve, useGMRES, verbose)
%Apply iterations to a multiscale problem
%
% SYNOPSIS:
%   [p_ms, report] = solveMultiscaleIteratively(A, q, basis, getSmootherFn, tol, iterations, LinSolve, useGMRES)
%
% DESCRIPTION:
%   Iteratively improve upon a solution using a multiscale preconditioner,
%   either as standalone or accelerated by GMRES.
%
% REQUIRED PARAMETERS:
%   A,q           - Linear system and right hand side. This solver attempts
%                   to solve Ap = q.
%
%   basis         - Basis functions as given by getMultiscaleBasis.
%
%   getSmootherFn - Smoother handle from getSmootherFunction.
%
%   tol           - Tolerance for convergence. A solution x is deemed to be
%                   converged if norm(A*x -q, 2)/norm(q) <= tol.
%
%   iterations    - Maximum number of iterations to apply.
%
%   LinSolve      - (OPTIONAL) Linear solver for coarse scale system.
%
%   useGMRES      - (OPTIONAL) Boolean indicating if GMRES should be used 
%                   to accelerate the solution process. Default false.
%
% RETURNS:
%   p_ms     - Solution, either converged or stopped by max iterations.
%
%   report   - Convergence report.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
    if nargin < 10
        verbose = mrstVerbose();
        if nargin < 9
            useGMRES = false;
            if nargin < 8
                LinSolve = @(A, b) mldivide(A, b);
            end
        end
    end
    B = basis.B;
    R = basis.R;
    
    A_c = R*A*B;    
    mssolver = @(d) B*LinSolve(A_c, (R*d));
    
    if isempty(p_ms)
        p_ms = mssolver(q);
    end
    
    if iterations > 0
        assert(~isempty(getSmootherFn), 'Need smoother function if iterations are to be used');
        
        % smoother is now fn handle of type @(defect)
        if nargin(getSmootherFn) == 1
            % We already have the smoother
            smoother = getSmootherFn;
        else
            % We have a setup phase
            smoother = getSmootherFn(A, q);
        end
        
        
        prec = @(b) twoStepMultiscalePreconditioner(A, b, mssolver, smoother);
        residuals = inf(iterations + 1, 1);
        
        if useGMRES
            [p_ms, flag, res, itcount, resvec] = gmres(A, q, [], tol, iterations, prec, [], p_ms);
            itcount = itcount(2);
            % GMRES gives out preconditioned residuals. We try to get the
            % actual values by setting the final preconditioned residual
            % equal to the final relative residual scaled by a constant.
            resvec = resvec.*(res./resvec(end));
        else
            [p_ms, flag, res, itcount, resvec] = simpleIterativeSolver(A, q, tol, iterations, prec, p_ms);
        end
        
        % Ensure residuals are always output as iterations + array, even if
        % the solver returned before it used up all the iterations. This
        % makes plotting etc much easier.
        residuals(1:numel(resvec)) = resvec;
        
        % Set lowest seen residual for the remaining of the vector
        [mr, x] = min(residuals);
        residuals(residuals > mr & (1:numel(residuals))' > x) = 0;

        if flag == 1 || flag == 2 || flag == 3
            itcount = iterations;
            dispif(verbose, 'Multiscale solver did not converge to desired precision\n');
        end
    else
        itcount = 0;
        res = norm(A*p_ms - q, 2)/norm(q, 2);
        residuals = res;
    end
    dispif(verbose, 'Final residual: %g after %d iterations (tolerance: %g)\n', res, itcount, tol);
    report.iterations = itcount;
    report.finalResidual = res;
    report.resvec = residuals;
    report.A_coarse = A_c;
end
