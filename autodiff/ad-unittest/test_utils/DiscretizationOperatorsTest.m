function tests = DiscretizationOperatorsTest
%Test suite for the diagonal backend's low-level discretization operators
%
% DESCRIPTION:
%   `discreteDivergence`, `singlePointUpwind`, `faceAverage` and
%   `twoPointGradient` (in `ad-core/backends/diagonal/operators/`) are
%   the flux/upwind/gradient kernels underlying nearly all flow physics
%   in every modern MRST model's `FlowDiscretization`, but had no
%   dedicated unit tests -- only indirect exercise via full-simulation
%   regression tests, always wrapped in AD. This suite checks their
%   plain-MATLAB (non-AD, non-mex) value computation against a hand
%   worked five-cell, four-interior-face 1D neighbor list, which is
%   exactly the code path used whenever MEX acceleration is
%   unavailable or disabled.
%
% SEE ALSO:
%   discreteDivergence, singlePointUpwind, faceAverage,
%   twoPointGradient, functiontests

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
   mrstModule reset
   mrstModule add ad-unittest ad-core

   % A 1D chain of 5 cells, 4 interior faces: face i connects cell i
   % and cell i+1.
   t.TestData.N  = [1 2; 2 3; 3 4; 4 5];
   t.TestData.v  = [10; 20; 30; 40; 50]; % per-cell values
   t.TestData.nc = 5;
end

%--------------------------------------------------------------------------

function testSinglePointUpwindSelectsFlagIndicatedNeighbor(t)    %#ok<DEFNU>
   % For each face, the upwind value must be v(N(:,1)) where flag is
   % true, and v(N(:,2)) where flag is false.
   N = t.TestData.N;
   v = t.TestData.v;

   flag = [true; false; true; false];
   vf = singlePointUpwind(flag, N, v, false);

   expected = [v(N(1,1)); v(N(2,2)); v(N(3,1)); v(N(4,2))];
   verifyEqual(t, vf, expected);
end

%--------------------------------------------------------------------------

function testFaceAverageIsArithmeticMean(t)                      %#ok<DEFNU>
   % Face-averaged value must be the arithmetic mean of the two
   % neighboring cell values.
   N = t.TestData.N;
   v = t.TestData.v;

   vf = faceAverage(N, v, false);

   expected = 0.5*(v(N(:,1)) + v(N(:,2)));
   verifyEqual(t, vf, expected);
   verifyEqual(t, vf, [15; 25; 35; 45]);
end

%--------------------------------------------------------------------------

function testTwoPointGradientIsDownwindMinusUpwind(t)            %#ok<DEFNU>
   % Discrete gradient across a face must be v(N(:,2)) - v(N(:,1)).
   N = t.TestData.N;
   v = t.TestData.v;

   g = twoPointGradient(N, v, [], false);

   expected = v(N(:,2)) - v(N(:,1));
   verifyEqual(t, g, expected);
   verifyEqual(t, g, [10; 10; 10; 10]);
end

%--------------------------------------------------------------------------

function testDiscreteDivergenceFluxOnly(t)                        %#ok<DEFNU>
   % With no accumulation term, the divergence at each cell must equal
   % the sum of fluxes leaving through faces where the cell is N(:,1)
   % minus the sum of fluxes entering through faces where the cell is
   % N(:,2) -- i.e. accumarray(N(:,1),flux) - accumarray(N(:,2),flux).
   N    = t.TestData.N;
   nc   = t.TestData.nc;
   flux = [10; 20; 30; 40]; % one value per interior face

   opts = struct('useMex', false, 'N', N, 'nc', nc);
   d = discreteDivergence([], flux, opts);

   expected = accumarray(N(:,1), flux, [nc, 1]) - accumarray(N(:,2), flux, [nc, 1]);
   verifyEqual(t, d, expected);
   verifyEqual(t, d, [10; 10; 10; 10; -40]);
end

%--------------------------------------------------------------------------

function testDiscreteDivergenceAddsAccumulationTerm(t)           %#ok<DEFNU>
   % With an accumulation term supplied, it must simply be added
   % on top of the flux divergence computed above.
   N    = t.TestData.N;
   nc   = t.TestData.nc;
   flux = [10; 20; 30; 40];
   acc  = [1; 2; 3; 4; 5];

   opts = struct('useMex', false, 'N', N, 'nc', nc);
   d = discreteDivergence(acc, flux, opts);

   expected = accumarray(N(:,1), flux, [nc, 1]) - accumarray(N(:,2), flux, [nc, 1]) + acc;
   verifyEqual(t, d, expected);
end
