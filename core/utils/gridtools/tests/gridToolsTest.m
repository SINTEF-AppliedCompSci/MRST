function tests = gridToolsTest
%Test suite for grid indexing helpers (gridCellNodes, gridCellFaces, gridCellNo)
%
% SEE ALSO:
%   gridCellNodes, gridCellFaces, gridCellNo, functiontests

%{
Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.

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

   tests = functiontests(localfunctions);
end

%--------------------------------------------------------------------------

function setupOnce(t)                                           %#ok<DEFNU>
   t.TestData.G = cartGrid([10, 5, 3]);
   t.TestData.c = [1, 2, 9, 7, 13, 51, 14, 40, 93];
end

%--------------------------------------------------------------------------

function testCellFaceSubsetIndexing(t)                          %#ok<DEFNU>
   % A subset of faces for a single cell, extracted through the
   % (cell, facemap) indirection array, must equal the faces obtained by
   % indexing directly into that cell's slice of the map.
   G = t.TestData.G;
   c = t.TestData.c;

   [faces, facemap] = gridCellFaces(G, c);

   target = 5;
   fi = facemap(target):facemap(target+1)-1;

   verifyEqual(t, faces(fi), indirectionSub(target, facemap, faces));
end

%--------------------------------------------------------------------------

function testCellNodesMatchesCoords(t)                          %#ok<DEFNU>
   % gridCellNodes must return exactly the node indices addressable
   % through G.nodes.coords for the requested cell subset.
   G = t.TestData.G;
   c = t.TestData.c;

   [nodes, nodemap] = gridCellNodes(G, c);
   verifyEqual(t, nodemap(end)-1, numel(nodes));
   verifyTrue(t, all(nodes >= 1 & nodes <= G.nodes.num));
end

%--------------------------------------------------------------------------

function testCellNoConsistency(t)                               %#ok<DEFNU>
   % Every half-face must belong to exactly one cell, and no cell may
   % appear in more than one contiguous block of the cellNo array.
   G = t.TestData.G;
   cellno = gridCellNo(G);

   fp = G.cells.facePos;
   seen = [];
   for i = 1:numel(fp)-1
      v = unique(cellno(fp(i):(fp(i+1)-1)));
      verifyEqual(t, numel(v), 1);
      verifyTrue(t, isempty(intersect(v, seen)));
      seen = [seen, v]; %#ok<AGROW>
   end
end
