function [x, flag, res, it, resvec] = simpleIterativeSolver(A, q, tol, it, prec, x)
    if nargin < 6
        x = zeros(size(q));
    end
    resvec = nan(it, 1);
    
    d = q - A*x;
    flag = 1;
    
    for it = 1:it
        res = norm(d, 2)/norm(q, 2);
        if res <= tol
            flag = 0;
            break
        end
        
        if it > 1 && resvec(it-1) < res
            flag = 3;
            warning('Convergence issues, aborting')
            break
        end
        x = x + prec(d);
        d = q - A*x;
        resvec(it) = res;
    end
    resvec = resvec(~isnan(resvec));
end
