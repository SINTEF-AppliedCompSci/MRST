function A = getIncomp1PhMatrix(G, T, state, fluid)
% Get TPFA-like incompressible system matrix via incompTPFA.

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
if nargin < 3
    state = initResSol(G, 0);
end

if nargin < 4
    require incomp
    fluid = initSingleFluid('rho', 1, 'mu', 1);
end

nf = size(getNeighbourship(G, 'topological', true), 1);
use_trans = numel(T) == nf;
state = incompTPFA(state, G, T, fluid, 'MatrixOutput', true, ...
                'LinSolve', @(A, x) 0*x, 'use_trans', use_trans);
A = state.A;
% Undo magic scaling from the inside of incompTPFA
A(1, 1) = A(1, 1)/2;
end
