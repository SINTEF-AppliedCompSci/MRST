function [p, solvetime] = reconstructPressureNormalized(CG, pressure, A, rhs)
% Solve reconstruction problem for multiscale methods

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
    D = formReconstructionMatrixLocal(A, CG.partition, CG.cells.centers);
    
%     tmp1 = warning('query','MATLAB:nearlySingularMatrix');
%     tmp2 = warning('query','MATLAB:singularMatrix');
%     warning('off','MATLAB:nearlySingularMatrix')
%     warning('off','MATLAB:singularMatrix')
    t = tic();
    p = mldivide(D, rhs - (A - D)*pressure);
    solvetime = toc(t);
%     warning(tmp1.state,'MATLAB:nearlySingularMatrix')
%     warning(tmp2.state,'MATLAB:singularMatrix')
end

function D = formReconstructionMatrixLocal(A, partition, centers)
    N = size(A, 1);
    [i, j, d] = find(A);
    isCenter = false(N, 1);
    isCenter(centers) = true;
    
    fix = isCenter(i) & isCenter(j);
    d(fix) = 2*d(fix);
    
    pi = partition(i);
    pj = partition(j);

    % Connections where both cells correspond to the same coarse block
    int = pi == pj;

    D  = sparse(i(int), j(int), d(int), N, N);
    D = D + diag(sum(A - D, 2));
end