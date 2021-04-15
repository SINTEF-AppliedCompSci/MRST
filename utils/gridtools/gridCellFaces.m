function [cf, p] = gridCellFaces(G, c)
% Find faces corresponding to a set of cells
%
% SYNOPSIS:
%   [cf, p] = gridCellFaces(G, c)
%
% PARAMETERS:
%   G    - Grid structure
%   c    - Cells where the fine faces are desired
%
%
% RETURNS:
%   p    - indirectionmap into `n`. The faces of cell `c(i)` is found at
%          positions `p(i):p(i+1)-1` in `n`
%   cf   - cell face positions in `G`

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


    % Find number of faces per cell
    nf = diff([G.cells.facePos(c), G.cells.facePos(c+1)], [],2);

    % Find the faces
    cf = G.cells.faces(mcolon(G.cells.facePos( c ), ...
                              G.cells.facePos(c+1) - 1),1);

    p = cumsum([1; double(nf)]);
end
