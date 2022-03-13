function [x, flag, res, itNo, resvec] = simpleIterativeSolver(A, q, tol, it, prec, x)
% Simple preconditioned iterative solver with same syntax as GMRES

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
    if nargin < 6
        x = zeros(size(q));
    end
    resvec = nan(it, 1);
    
    d = q - A*x;
    flag = 1;
    
    for itNo = 1:it
        res = norm(d, 2)/norm(q, 2);
        if res <= tol
            flag = 0;
            break
        end
        
        if it > 1 && 2*resvec(it-1) < res
            flag = 3;
            warning('Convergence issues, aborting')
            break
        end
        x = x + prec(d);
        d = q - A*x;
        resvec(itNo) = res;
    end
    resvec = resvec(~isnan(resvec));
end
