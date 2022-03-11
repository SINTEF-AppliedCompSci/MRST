function ix = searchForBoundaryFaces(G, direction, i1, i2, i3)
%Search for boundary faces within (a subset of) a logically Cartesian grid
%
% SYNOPSIS:
%   ix = searchForBoundaryFaces(G, side)
%   ix = searchForBoundaryFaces(G, side, i1, i2, i3)
%
% PARAMETERS:
%   G     - Grid data structure.
%
%   side  - Global side from which to extract face indices.  String.  Must
%           (case insensitively) match one of six alias groups:
%
%              1) {'West' , 'XMin', 'Left'  }
%              2) {'East' , 'XMax', 'Right' }
%              3) {'South', 'YMin', 'Back'  }
%              4) {'North', 'YMax', 'Front' }
%              5) {'Upper', 'ZMin', 'Top'   }
%              6) {'Lower', 'ZMax', 'Bottom'}
%
%           These groups correspond to the cardinal directions mentioned as
%           the first alternative in each group.
%
%   i1,i2,i3 - Index ranges in which one is to look for boundary faces in
%           the three axial directions
%
% RETURNS:
%   ix    - Required face indices.

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

% $Date: 2012-01-30 11:41:03 +0100 (Mon, 30 Jan 2012) $
% $Revision: 9020 $

if ~exist('i1', 'var'), i1 = []; end
if ~exist('i2', 'var'), i2 = []; end
if ~exist('i3', 'var'), i3 = []; end

i1 = reshape(i1, [], 1);
i2 = reshape(i2, [], 1);
i3 = reshape(i3, [], 1);

switch lower(direction)
   case {'left'  , 'west' , 'xmin'}
      faceTag = 1;
   case {'right' , 'east' , 'xmax'}
      faceTag = 2;
   case {'back'  , 'south', 'ymin'}
      faceTag = 3;
   case {'front' , 'north', 'ymax'}
      faceTag = 4;
   case {'top'   , 'upper', 'zmin'}
      faceTag = 5;
   case {'bottom', 'lower', 'zmax'}
      faceTag = 6;
   otherwise
      error(1,'Illigal direction');
end
[Ix{1:3}] = ind2sub(G.cartDims,G.cells.indexMap);
cellno = rldecode(1:G.cells.num, diff(G.cells.facePos),2)';
i      = any(G.faces.neighbors ==0, 2);
j      = false(max(G.cells.faces(:,2)),1);
if isempty(i1)
   pick_I = true(G.cartDims(1),1);
else
   pick_I = false(G.cartDims(1),1); pick_I(i1) = true;
end
if isempty(i2)
   pick_J = true(G.cartDims(2),1);
else
   pick_J = false(G.cartDims(2),1); pick_J(i2) = true;
end
if numel(G.cartDims) == 2
    pick_K = true;
else
    assert(numel(G.cartDims) == 3);
    if isempty(i3) && (numel(G.cartDims) == 3)
        pick_K = true(G.cartDims(3),1);
    else
        pick_K = false(G.cartDims(3),1); pick_K(i3) = true;
    end
end

j(faceTag) = true;
p  = pick_I(Ix{1}) & pick_J(Ix{2}) & pick_K(Ix{3});
I  = i(G.cells.faces(:,1)) & j(G.cells.faces(:,2)) & p(cellno);
ix = G.cells.faces(I,1);
end
