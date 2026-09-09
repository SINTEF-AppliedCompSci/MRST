function tests = AddWellTest
%Test suite for addWell (well index, geometry and structural defaults)
%
% DESCRIPTION:
%   `addWell` (280+ references repo-wide) had no dedicated unit tests.
%   This suite checks its computed Peaceman well index against the
%   closed-form formula for a vertical well through a Cartesian cell,
%   both for isotropic and areally-anisotropic permeability, plus a
%   handful of structural defaults (naming, sign, status, depth offset)
%   that many downstream well-handling routines rely on.
%
% SEE ALSO:
%   addWell, computeWellIndex, functiontests

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
   % A 3x3x2 grid with dx=dy=10, dz=5, so the Peaceman drainage radius
   % and well index for a vertical ('z') well are easy to check by hand.
   G = cartGrid([3, 3, 2], [30, 30, 10]);
   t.TestData.G  = computeGeometry(G);
   t.TestData.dx = 10;
   t.TestData.dy = 10;
   t.TestData.dz = 5;
   t.TestData.rw = 0.1;
end

%--------------------------------------------------------------------------

function testPeacemanWellIndexIsotropic(t)                       %#ok<DEFNU>
   % For isotropic permeability k, a vertical well's index must match
   % the standard Peaceman formula
   %   re = 0.14*sqrt(dx^2 + dy^2),  WI = 2*pi*k*dz / log(re/rw).
   G  = t.TestData.G;
   dx = t.TestData.dx; dy = t.TestData.dy; dz = t.TestData.dz;
   rw = t.TestData.rw;
   k  = 100*milli*darcy;

   rock = struct('perm', repmat(k, G.cells.num, 3));
   W = addWell([], G, rock, 5, 'Radius', rw, 'Dir', 'z');

   re = 0.14*sqrt(dx^2 + dy^2);
   expectedWI = 2*pi*k*dz / log(re/rw);

   verifyEqual(t, W.WI, expectedWI, 'RelTol', 1e-10);
end

%--------------------------------------------------------------------------

function testPeacemanWellIndexAnisotropic(t)                     %#ok<DEFNU>
   % For areally-anisotropic permeability (kx ~= ky), the well index
   % must match the generalized Peaceman formula using the geometric
   % mean permeability ke = sqrt(kx*ky) and the anisotropic drainage
   % radius.
   G  = t.TestData.G;
   dx = t.TestData.dx; dy = t.TestData.dy; dz = t.TestData.dz;
   rw = t.TestData.rw;
   kx = 200*milli*darcy;
   ky = 50*milli*darcy;
   kz = 100*milli*darcy;

   rock = struct('perm', repmat([kx, ky, kz], G.cells.num, 1));
   W = addWell([], G, rock, 5, 'Radius', rw, 'Dir', 'z');

   wc  = 0.14;
   k21 = ky/kx; k12 = kx/ky;
   re  = 2*wc*sqrt(dx^2*sqrt(k21) + dy^2*sqrt(k12)) / (nthroot(k21,4) + nthroot(k12,4));
   ke  = sqrt(kx*ky);
   expectedWI = 2*pi*ke*dz / log(re/rw);

   verifyEqual(t, W.WI, expectedWI, 'RelTol', 1e-10);
end

%--------------------------------------------------------------------------

function testDepthOffsetMatchesCentroidDifference(t)             %#ok<DEFNU>
   % `dZ` must be each perforated cell's centroid depth measured from
   % the depth of the shallowest perforated cell's TOP contact (its
   % centroid minus half its own thickness dz), per addWell's
   % documented definition of the reference ("highest horizontal
   % contact") depth.
   G = t.TestData.G;
   rock = struct('perm', repmat(100*milli*darcy, G.cells.num, 3));

   cells = [5, 14]; % same (i,j), two different k-layers
   W = addWell([], G, rock, cells, 'Radius', 0.1, 'Dir', 'z');

   z      = G.cells.centroids(cells, 3);
   topZ   = z(1) - t.TestData.dz/2;
   verifyEqual(t, W.dZ, z - topZ, 'AbsTol', 1e-10);
end

%--------------------------------------------------------------------------

function testStructuralDefaultsAndMultipleWells(t)                %#ok<DEFNU>
   % Adding a second well must append (not replace) the well array,
   % auto-name wells sequentially, default open status/completions and
   % the water-injection composition, and derive an undetermined sign
   % (0) for a bhp-controlled well but a definite sign from a
   % rate-controlled well's target value.
   G = t.TestData.G;
   rock = struct('perm', repmat(100*milli*darcy, G.cells.num, 3));

   W = addWell([], G, rock, 5, 'Type', 'bhp', 'Val', 100*barsa, ...
               'Radius', 0.1, 'Dir', 'z');
   W = addWell(W, G, rock, 14, 'Type', 'rate', 'Val', 10/day, ...
               'Radius', 0.1, 'Dir', 'z');

   verifyEqual(t, numel(W), 2);
   verifyEqual(t, W(1).name, 'W1');
   verifyEqual(t, W(2).name, 'W2');
   verifyEqual(t, W(1).sign, 0);  % bhp control: sign undetermined
   verifyEqual(t, W(2).sign, 1);  % positive rate => injector
   verifyEqual(t, W(1).status, true);
   verifyEqual(t, W(1).cstatus, true);
   verifyEqual(t, W(1).compi, [1, 0, 0]);
end
