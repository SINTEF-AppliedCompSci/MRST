function tests = NonLinearSolverHelpersTest
%Test suite for NonLinearSolver's residual-history diagnostics
%
% DESCRIPTION:
%   `NonLinearSolver` drives the Newton loop used by every model, but
%   is otherwise only exercised end-to-end via full-simulation
%   regression tests. `checkForOscillations` and `checkForStagnation`
%   are two of its pure, model-independent helpers (given only a
%   residual-history matrix and an iteration index), used to trigger
%   line search / stabilization -- yet had no dedicated unit test.
%
% SEE ALSO:
%   NonLinearSolver, functiontests

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
   t.TestData.solver = NonLinearSolver();
end

%--------------------------------------------------------------------------

function testOscillationDetectedForAlternatingResiduals(t)       %#ok<DEFNU>
   % A perfectly alternating residual (1,2,1,2,...) must be flagged as
   % oscillating from the third entry onwards (the check needs three
   % points to form a forward/backward difference ratio), and never
   % before that.
   solver = t.TestData.solver;
   res = [1; 2; 1; 2];

   verifyEqual(t, solver.checkForOscillations(res, 1), false);
   verifyEqual(t, solver.checkForOscillations(res, 2), false);
   verifyEqual(t, solver.checkForOscillations(res, 3), true);
   verifyEqual(t, solver.checkForOscillations(res, 4), true);
end

%--------------------------------------------------------------------------

function testNoOscillationForMonotonicallyDecreasingResiduals(t) %#ok<DEFNU>
   % A monotonically (geometrically) decreasing residual must never be
   % flagged as oscillating.
   solver = t.TestData.solver;
   res = [100; 50; 25; 12.5; 6.25];

   for idx = 3:numel(res)
      verifyEqual(t, solver.checkForOscillations(res, idx), false);
   end
end

%--------------------------------------------------------------------------

function testStagnationDetectedForSlowlyChangingResiduals(t)     %#ok<DEFNU>
   % A residual that changes by less than the (default) stagnation
   % tolerance between consecutive iterations must be flagged as
   % stagnated; one that is still dropping quickly must not.
   solver = t.TestData.solver;
   verifyEqual(t, solver.stagnateTol, 0.01, 'AbsTol', 1e-12);

   res = [10; 5; 5.001; 5.0015];

   verifyEqual(t, solver.checkForStagnation(res, 1), false);
   verifyEqual(t, solver.checkForStagnation(res, 2), false); % 5 vs 10: big relative change
   verifyEqual(t, solver.checkForStagnation(res, 3), true);  % 5.001 vs 5: tiny relative change
   verifyEqual(t, solver.checkForStagnation(res, 4), true);  % 5.0015 vs 5.001: tiny relative change
end

%--------------------------------------------------------------------------

function testDiagnosticsOperateColumnwiseOnMultipleResiduals(t)  %#ok<DEFNU>
   % With more than one residual tracked side by side (one per
   % column), both diagnostics must operate independently per column.
   solver = t.TestData.solver;
   % Column 1 oscillates, column 2 decreases monotonically.
   res = [1, 100; 2, 50; 1, 25];

   verifyEqual(t, solver.checkForOscillations(res, 3), [true, false]);
end
