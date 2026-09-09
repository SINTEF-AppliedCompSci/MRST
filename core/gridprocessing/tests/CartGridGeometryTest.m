function tests = CartGridGeometryTest
%Test suite for cartGrid/tensorGrid geometry via computeGeometry
%
% DESCRIPTION:
%   `cartGrid`, `tensorGrid` and `computeGeometry` are, by reference
%   count, the most heavily depended-on functions in all of MRST, yet
%   have no dedicated unit tests: they are only exercised indirectly
%   through whatever example or regression test happens to build a
%   grid. This suite checks their output against closed-form geometry
%   for a small, non-cubic Cartesian grid where cell volumes, face
%   areas/normals and centroids are all easy to compute by hand.
%
% SEE ALSO:
%   cartGrid, tensorGrid, computeGeometry, functiontests

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

function setupOnce(t)                                            %#ok<DEFNU>
   % A 2-by-3-by-1 grid with anisotropic cell size dx=2, dy=3, dz=2, so
   % that cell volumes, per-direction face areas and centroids are all
   % distinct and easy to verify by hand.
   t.TestData.celldim = [2, 3, 1];
   t.TestData.physdim = [4, 9, 2];
   t.TestData.dx      = [2, 3, 2];
   G = cartGrid(t.TestData.celldim, t.TestData.physdim);
   t.TestData.G = computeGeometry(G);
end

%--------------------------------------------------------------------------

function testCellDimsAndCount(t)                                 %#ok<DEFNU>
   % The grid must report the requested logical dimensions and the
   % corresponding number of cells.
   G = t.TestData.G;
   verifyEqual(t, G.cartDims, t.TestData.celldim);
   verifyEqual(t, G.cells.num, prod(t.TestData.celldim));
end

%--------------------------------------------------------------------------

function testCellVolumes(t)                                      %#ok<DEFNU>
   % Every cell has the same volume dx*dy*dz, and the volumes must sum
   % to the total physical volume of the domain.
   G  = t.TestData.G;
   dx = t.TestData.dx;

   expectedVolume = prod(dx);
   verifyEqual(t, G.cells.volumes, repmat(expectedVolume, G.cells.num, 1), ...
      'AbsTol', 1e-10);
   verifyEqual(t, sum(G.cells.volumes), prod(t.TestData.physdim), 'AbsTol', 1e-10);
end

%--------------------------------------------------------------------------

function testCellCentroids(t)                                    %#ok<DEFNU>
   % Cell centroids must sit at the midpoint of each logical (i,j,k)
   % cell, i.e. ((i-1/2)*dx, (j-1/2)*dy, (k-1/2)*dz).
   G  = t.TestData.G;
   dx = t.TestData.dx;

   [ii, jj, kk] = ind2sub(G.cartDims, (1:G.cells.num)');
   expected = [(ii - 0.5)*dx(1), (jj - 0.5)*dx(2), (kk - 0.5)*dx(3)];

   verifyEqual(t, G.cells.centroids, expected, 'AbsTol', 1e-10);
end

%--------------------------------------------------------------------------

function testFaceAreasByDirection(t)                              %#ok<DEFNU>
   % Face area must equal the product of the two cell dimensions
   % orthogonal to the face's normal direction: dy*dz for x-facing
   % faces, dx*dz for y-facing, dx*dy for z-facing.
   G  = t.TestData.G;
   dx = t.TestData.dx;

   n = G.faces.normals;
   isXFace = n(:,1) ~= 0 & n(:,2) == 0 & n(:,3) == 0;
   isYFace = n(:,2) ~= 0 & n(:,1) == 0 & n(:,3) == 0;
   isZFace = n(:,3) ~= 0 & n(:,1) == 0 & n(:,2) == 0;

   % Every face must be axis-aligned in exactly one direction.
   verifyTrue(t, all(isXFace | isYFace | isZFace));
   verifyEqual(t, nnz(isXFace & isYFace) + nnz(isYFace & isZFace) + ...
                   nnz(isXFace & isZFace), 0);

   verifyEqual(t, G.faces.areas(isXFace), ...
      repmat(dx(2)*dx(3), nnz(isXFace), 1), 'AbsTol', 1e-10);
   verifyEqual(t, G.faces.areas(isYFace), ...
      repmat(dx(1)*dx(3), nnz(isYFace), 1), 'AbsTol', 1e-10);
   verifyEqual(t, G.faces.areas(isZFace), ...
      repmat(dx(1)*dx(2), nnz(isZFace), 1), 'AbsTol', 1e-10);
end

%--------------------------------------------------------------------------

function testFaceNormalMagnitudeEqualsArea(t)                     %#ok<DEFNU>
   % computeTrans relies on face normals having magnitude equal to the
   % face area (documented invariant of computeGeometry).
   G = t.TestData.G;
   normalMagnitude = sqrt(sum(G.faces.normals.^2, 2));
   verifyEqual(t, normalMagnitude, G.faces.areas, 'AbsTol', 1e-10);
end
