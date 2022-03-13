function D = formReconstructionMatrix(A, partition, keepDiag)
% Form matrix for multiscale flux reconstruction

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
    if nargin == 2
        keepDiag = false;
    end
    N = size(A, 1);
    [i, j, d] = find(A);

    pi = partition(i);
    pj = partition(j);

    % Connections where both cells correspond to the same coarse block
    int = pi == pj;

    D  = sparse(i(int), j(int), d(int), N, N);
    if ~keepDiag
        D = D + diag(sum(A - D, 2));
    end
end
