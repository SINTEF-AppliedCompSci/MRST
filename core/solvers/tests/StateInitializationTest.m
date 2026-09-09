function tests = StateInitializationTest
%Test suite for initResSol, initWellSol and initState
%
% DESCRIPTION:
%   `initResSol` (326 references), `initState` (181) and `initWellSol`
%   (75) set up the reservoir/well state structure nearly every solver
%   and example builds on, but had no dedicated unit tests. This suite
%   checks their broadcasting rules, output shapes, and the one
%   consistency check `initState` performs between well compositions
%   and reservoir phase count.
%
% SEE ALSO:
%   initResSol, initWellSol, initState, addWell, functiontests

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
   G = cartGrid([3, 3, 2], [30, 30, 10]);
   t.TestData.G    = computeGeometry(G);
   t.TestData.rock = struct('perm', repmat(100*milli*darcy, t.TestData.G.cells.num, 3));
end

%--------------------------------------------------------------------------

function testInitResSolBroadcastsScalarsAndDefaultsSinglePhase(t) %#ok<DEFNU>
   % A scalar pressure must be broadcast to every cell, flux must be
   % all-zero with one entry per face, and the default saturation (no
   % s0 given) must be full (single-phase) saturation in every cell.
   G = t.TestData.G;
   p0 = 200*barsa;

   state = initResSol(G, p0);

   verifyEqual(t, state.pressure, repmat(p0, G.cells.num, 1));
   verifyEqual(t, state.flux, zeros(G.faces.num, 1));
   verifyEqual(t, state.s, ones(G.cells.num, 1));
end

%--------------------------------------------------------------------------

function testInitResSolBroadcastsMultiphaseSaturation(t)         %#ok<DEFNU>
   % A 1-by-np saturation vector must be broadcast to every cell, and a
   % per-cell pressure vector of the right size must pass through
   % unchanged.
   G = t.TestData.G;
   p0 = (100:100+G.cells.num-1)'*barsa;
   s0 = [0.2, 0.3, 0.5];

   state = initResSol(G, p0, s0);

   verifyEqual(t, state.pressure, p0);
   verifyEqual(t, state.s, repmat(s0, G.cells.num, 1));
end

%--------------------------------------------------------------------------

function testInitResSolRejectsMismatchedSaturationRows(t)        %#ok<DEFNU>
   % A saturation array whose row count is neither 1 nor G.cells.num
   % must be rejected.
   G = t.TestData.G;
   verifyError(t, @() initResSol(G, 100*barsa, ones(3, 2)), ...
      'initResSol:InitSat:Inconsistent');
end

%--------------------------------------------------------------------------

function testInitWellSolPerWellShapes(t)                          %#ok<DEFNU>
   % Each well's flux must be all-zero with one entry per perforation,
   % and its pressure must equal the supplied uniform initial pressure.
   G = t.TestData.G;
   rock = t.TestData.rock;
   p0 = 150*barsa;

   W = addWell([], G, rock, 5, 'Radius', 0.1, 'Dir', 'z');
   W = addWell(W, G, rock, [8, 14], 'Radius', 0.1, 'Dir', 'z');

   wellSol = initWellSol(W, p0);

   verifyEqual(t, numel(wellSol), 2);
   verifyEqual(t, wellSol(1).flux, zeros(1, 1));
   verifyEqual(t, wellSol(2).flux, zeros(2, 1));
   verifyEqual(t, wellSol(1).pressure, p0);
   verifyEqual(t, wellSol(2).pressure, p0);
end

%--------------------------------------------------------------------------

function testInitStateCombinesReservoirAndWellState(t)           %#ok<DEFNU>
   % initState must produce both the reservoir fields (as initResSol
   % would) and a wellSol field (as initWellSol would) when wells are
   % supplied, using the same initial pressure for both.
   G = t.TestData.G;
   rock = t.TestData.rock;
   p0 = 175*barsa;
   s0 = [0.2, 0.3, 0.5];

   W = addWell([], G, rock, 5, 'Radius', 0.1, 'Dir', 'z', 'compi', s0);

   state = initState(G, W, p0, s0);

   verifyEqual(t, state.pressure, repmat(p0, G.cells.num, 1));
   verifyEqual(t, state.s, repmat(s0, G.cells.num, 1));
   verifyEqual(t, state.wellSol.pressure, p0);
   verifyEqual(t, state.wellSol.flux, zeros(1, 1));
end

%--------------------------------------------------------------------------

function testInitStateRejectsWellReservoirPhaseMismatch(t)       %#ok<DEFNU>
   % initState must reject an s0 whose phase count does not match the
   % number of columns in the wells' injection composition (compi).
   G = t.TestData.G;
   rock = t.TestData.rock;

   W = addWell([], G, rock, 5, 'Radius', 0.1, 'Dir', 'z', ...
               'compi', [1, 0, 0]); % 3-phase composition

   verifyError(t, @() initState(G, W, 100*barsa, [0.2, 0.8]), ...
      ?MException); % 2-phase s0
end
