function b = computeAdjointRHS(G, W, f_res)
% Compute adjoint 'pressure' rhs
% for use in e.g. function solveAdjointPressureSystem
%
% SYNOPSIS:
%   b = computeAdjoint(G, W, f_res, f_w)
%
% DESCRIPTION:
% Computes adjoint 'pressure' rhs as input to solveIncompFlow
% in function solveAdjointPressureSystem
%
% PARAMETERS:
%   G     - Grid data structure.
%
%   W     - Well structure as defined by addWell &c.
%
%   f_res - Adjoint reservoir 'pressure' condtions
%
% RETURNS:
%   b     - Ajoint pressure rhs to be passed directly as option
%           'rhs' to solveIncompFlow.

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


assert( numel(f_res) == size(G.cells.faces,1));
f_w = [];

% unpack adjoint well 'pressure' conditions
if ~isempty(W),
   S   = [ W.S   ];
   RHS = [ S.RHS ];
   f_w = {RHS.f};
   f_w = vertcat(f_w{:});
   h_w = {RHS.h};
   h_w = vertcat(h_w{:});
end

% b = [f_res; f_w; g; h_res; h_w]
b    = cell([1, 3]);
b{1} = vertcat(f_res, f_w);
b{2} =         zeros([G.cells.num, 1]);
b{3} = vertcat(zeros([G.faces.num, 1]), h_w);
end
