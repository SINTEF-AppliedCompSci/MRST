function [p, solvetime] = reconstructPressure(partition, pressure, A, rhs)
% Solve reconstruction problem for multiscale methods

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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
    D = formReconstructionMatrix(A, partition);
    n = size(A, 1);
    % Add a tiny bit of noise to overcome a Matlab UMFPACK issue in MRST
    % 2018b and onwards
    d = diag(A);
    D = D - sparse(1:n, 1:n, min(d)*1e-8, n, n);
    tmp1 = warning('query','MATLAB:nearlySingularMatrix');
    tmp2 = warning('query','MATLAB:singularMatrix');
    warning('off','MATLAB:nearlySingularMatrix')
    warning('off','MATLAB:singularMatrix')
    t = tic();
    p = mldivide(D, rhs - (A - D)*pressure);
    assert(all(isfinite(p)))
    solvetime = toc(t);
    warning(tmp1.state,'MATLAB:nearlySingularMatrix')
    warning(tmp2.state,'MATLAB:singularMatrix')
end
