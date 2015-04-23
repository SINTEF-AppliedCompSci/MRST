function [p_ms, report] = solveMultiscaleIteratively(A, q, basis, getSmootherFn, tol, iterations, LinSolve, useGMRES)
    if nargin < 8
        useGMRES = false;
        if nargin < 7
            LinSolve = @(A, b) mldivide(A, b);
        end
    end
    B = basis.B;
    R = basis.R;
    
    A_c = R*A*B;    
    mssolver = @(d) B*LinSolve(A_c, (R*d));
    
    p_ms = mssolver(q);
    
    if iterations > 0
        assert(~isempty(getSmootherFn), 'Need smoother function if iterations are to be used');
        
        % smoother is now fn handle of type @(defect)
        smoother = getSmootherFn(A, q);
        
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
            dispif(mrstVerbose(), 'Multiscale solver did not converge to desired precision\n');
        end
    else
        itcount = 0;
        res = norm(A*p_ms - q, 2)/norm(q, 2);
        residuals = res;
    end
    dispif(mrstVerbose, 'Converged to %g after %d iterations\n', res, itcount);
    report.iterations = itcount;
    report.finalResidual = res;
    report.resvec = residuals;
    report.A_coarse = A_c;
end
