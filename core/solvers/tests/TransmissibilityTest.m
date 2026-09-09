function tests = TransmissibilityTest
%Test suite for computeTrans and getFaceTransmissibility
%
% DESCRIPTION:
%   `computeTrans` (half-transmissibilities) and
%   `getFaceTransmissibility` (per-face transmissibilities, combining
%   the two one-sided contributions harmonically) are two of the most
%   heavily depended-on functions in `core/`, but had no dedicated unit
%   tests. This suite checks both against closed-form values for a
%   minimal two-cell Cartesian grid with anisotropic permeability,
%   where the half-transmissibility in each direction is
%   `k_dir * A_dir / (0.5*dx_dir)`.
%
% SEE ALSO:
%   computeTrans, getFaceTransmissibility, computeGeometry, functiontests

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
   % Two cells side by side in x, dx=dy=dz=2, with anisotropic
   % permeability [kx, ky, kz] = [2, 4, 6] so every direction gives a
   % distinct, easily hand-checked half-transmissibility:
   %   half-trans_x = kx*(dy*dz)/(dx/2) = 2*4/1  =  8
   %   half-trans_y = ky*(dx*dz)/(dy/2) = 4*4/1  = 16
   %   half-trans_z = kz*(dx*dy)/(dz/2) = 6*4/1  = 24
   G = cartGrid([2, 1, 1], [4, 2, 2]);
   t.TestData.G    = computeGeometry(G);
   t.TestData.rock = struct('perm', repmat([2, 4, 6], G.cells.num, 1));
   t.TestData.halfTransByDirection = [8, 8, 16, 16, 24, 24]; % codes 1..6
end

%--------------------------------------------------------------------------

function testHalfTransmissibilityByDirection(t)                  %#ok<DEFNU>
   % computeTrans must return, for every (cell, local face), the
   % one-sided transmissibility for that face's direction code
   % (1=W,2=E,3=S,4=N,5=T,6=B), independent of whether the face is
   % an interior or a boundary face.
   G    = t.TestData.G;
   rock = t.TestData.rock;

   T = computeTrans(G, rock);
   expected = t.TestData.halfTransByDirection(G.cells.faces(:,2))';

   verifyEqual(t, T, expected, 'AbsTol', 1e-10);
end

%--------------------------------------------------------------------------

function testFaceTransmissibilityHarmonicAverage(t)              %#ok<DEFNU>
   % getFaceTransmissibility must combine the two one-sided
   % transmissibilities of an interior face harmonically
   % (1/T = 1/T1 + 1/T2); for a boundary face (only one neighbor) it
   % must reduce to that single one-sided value.
   G    = t.TestData.G;
   rock = t.TestData.rock;

   Tf = getFaceTransmissibility(G, rock);
   n  = G.faces.normals;

   isInterior = all(G.faces.neighbors > 0, 2);
   isXFace = n(:,1) ~= 0 & n(:,2) == 0 & n(:,3) == 0;
   isYFace = n(:,2) ~= 0 & n(:,1) == 0 & n(:,3) == 0;
   isZFace = n(:,3) ~= 0 & n(:,1) == 0 & n(:,2) == 0;

   % The only interior face is the single x-face between the two cells;
   % its two one-sided contributions are both 8, giving a harmonic
   % average of 8*8/(8+8) = 4.
   verifyEqual(t, nnz(isInterior), 1);
   verifyTrue(t, isXFace(isInterior));
   verifyEqual(t, Tf(isInterior), 4, 'AbsTol', 1e-10);

   % Every boundary face reduces to its single one-sided transmissibility.
   isBoundary = ~isInterior;
   verifyEqual(t, Tf(isBoundary & isXFace), ...
      repmat(8,  nnz(isBoundary & isXFace), 1), 'AbsTol', 1e-10);
   verifyEqual(t, Tf(isBoundary & isYFace), ...
      repmat(16, nnz(isBoundary & isYFace), 1), 'AbsTol', 1e-10);
   verifyEqual(t, Tf(isBoundary & isZFace), ...
      repmat(24, nnz(isBoundary & isZFace), 1), 'AbsTol', 1e-10);
end

%--------------------------------------------------------------------------

function testScalarPermeabilityIsIsotropicShorthand(t)           %#ok<DEFNU>
   % A single-column (scalar) permeability must be equivalent to an
   % isotropic diagonal tensor with that value repeated in every
   % direction.
   G = t.TestData.G;
   k = 3;
   rockScalar = struct('perm', repmat(k, G.cells.num, 1));
   rockDiag   = struct('perm', repmat([k, k, k], G.cells.num, 1));

   Tscalar = computeTrans(G, rockScalar);
   Tdiag   = computeTrans(G, rockDiag);

   verifyEqual(t, Tscalar, Tdiag, 'AbsTol', 1e-10);
end
