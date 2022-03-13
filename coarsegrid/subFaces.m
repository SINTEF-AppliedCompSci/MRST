function [nsub, sub] = subFaces(g, cg)
%Extract fine-grid faces constituting individual coarse grid faces.
%
% SYNOPSIS:
%   [nsub, sub] = subFaces(G, CG)
%
% PARAMETERS:
%   G  - Grid data structure as described by grid_structure.
%   CG - Coarse grid data structure.
%
% RETURNS:
%   nsub - Number of fine-grid faces belonging to each individual coarse
%          grid face.  Specifically, nsub(i) is the number of fine-grid
%          faces belonging to coarse face 'i'.
%
%   sub  - Actual fine-grid subfaces represented in a condensed storage
%          format.  Specifically, if IX = CUMSUM([0; nsub]), then the
%          subfaces of coarse face 'i' are sub(IX(i) + 1 : IX(i + 1)).
%
% EXAMPLE:
%   % Generate moderately large grid and corresponding coarse grid
%   G  = cartGrid([60, 220, 85]);
%   p  = partitionUI(G, [5, 11, 17]);
%   CG = generateCoarseGrid(G, p);
%
%   % Extract subfaces
%   [nsub, sub] = subFaces(G, CG);
%
%   % Find fine-faces belonging to coarse face 10
%   sub_ix = cumsum([0; nsub]);
%   sf10   = sub(sub_ix(10) + 1 : sub_ix(10 + 1))
%
%   % Find index of coarse face to which fine-face 7 belongs
%   f2c = sparse(sub, 1, rldecode((1 : CG.faces.num) .', nsub));
%   cf7 = f2c(7)
%
% SEE ALSO:
%   `grid_structure`, `generateCoarseGrid`.

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


   assert (isfield(cg.faces, 'fconn'), ...
          ['This edition of function ''%s'' is only supported when ', ...
           'accompanied by the revised edition of function ', ...
           '''generateCoarseGrid'''], mfilename);

   assert (max(cg.faces.fconn) <= g.faces.num, ...
          ['The coarse grid does not appear to match the ', ...
           'fine-scale grid.']);

   nsub = diff(cg.faces.connPos);
   sub  = cg.faces.fconn;
end
