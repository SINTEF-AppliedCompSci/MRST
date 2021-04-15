function G_new = flipGrid(G)
% Flip a grid (z->x, x->y, y->z)
%
% SYNOPSIS:
%   function G_new = flipGrid(G);
%
% DESCRIPTION: Flip the grid (z->x, x->y, y->z)
%
% PARAMETERS:
%   G - Grid structure
%
% RETURNS:
%   G_new - Grid structure afte flipping

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

    G_new = G;
    per = [3, 1, 2];
    tags = reshape(1:6, 2, []);
    new_tag = reshape(tags(:, [2, 3, 1]), [], 1);
    G_new.nodes.coords = G_new.nodes.coords(:, per);
    G_new.cells.faces(:, 2) = new_tag(G.cells.faces(:, 2));
end
