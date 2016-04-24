function [p_ms, report] = blockSolverMultiscaleCycle(state, CG, fluid, A, q, p_ms, msrsb, tol, iterations, W, nearWell, verbose)
G = CG.parent;
nw = numel([W.cells]);
x = p_ms;
if iterations > 0
    resvec = nan(iterations, 1);
    d = q - A*x;
    flag = 1;
    residuals = inf(iterations + 1, 1);
    for i = 1:iterations
        res = norm(d, 2)/norm(q, 2);
        if res <= tol
            flag = 0;
            break
        end
        if i > 1 && (resvec(i-1) - res) < -eps*1000
            flag = 3;
            warning('Convergence issues, aborting')
            break
        end
        %-------------------------------------------------------------%
        [Ab, bb] = extractSubBlock(CG, state, fluid, x, nearWell, W);
        wellb = size(x)-nw+1:size(x);
        block = [nearWell;wellb'];
        x(block) = x(block) + Ab\(bb - Ab*x(block));
        x = x + msrsb(q - A*x);
        %-------------------------------------------------------------%
        d = q - A*x;
        resvec(i) = res;
    end
    p_ms = x;
    itcount = i;
    resvec = resvec(~isnan(resvec));
    
    % Ensure residuals are always output as iterations + array, even if
    % the solver returned before it used up all the iterations. This
    % makes plotting etc much easier.
    residuals(1:numel(resvec)) = resvec;
    
    % Set lowest seen residual for the remaining of the vector
    [mr, x] = min(residuals);
    residuals(residuals > mr & (1:numel(residuals)).' > x) = 0;
    
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
