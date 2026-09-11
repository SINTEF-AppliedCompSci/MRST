function tests = DiagonalBackendConsistencyTest
%Test suite checking DiagonalAutoDiffBackend agrees with SparseAutoDiffBackend
%
% DESCRIPTION:
%   `DiagonalAutoDiffBackend` is the default, high-performance AD backend
%   used by essentially every modern MRST model, but its diagonal
%   Jacobian representation is never compared directly against the
%   simple, easy-to-trust `SparseAutoDiffBackend` reference
%   implementation at the unit level -- only indirectly, through
%   full-simulation regression tests. This suite builds the same
%   expressions with both backends and checks that values and (once
%   converted to plain sparse matrices) Jacobians agree.
%
% SEE ALSO:
%   DiagonalAutoDiffBackend, SparseAutoDiffBackend, GenericAD, functiontests

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

   rng(20260909);
   t.TestData.xv = rand(6, 1) + 0.2;
   t.TestData.yv = rand(6, 1) + 0.2;
end

%--------------------------------------------------------------------------

function testArithmeticExpressionMatchesSparseBackend(t)         %#ok<DEFNU>
   % A composite expression exercising +, -, .*, ./, .^, exp and sin
   % together must give identical values and (fully assembled) Jacobians
   % under both backends.
   xv = t.TestData.xv;
   yv = t.TestData.yv;

   [xs, ys] = getBackendVars('sparse', xv, yv);
   zs = 2*xs - xs.*ys + xs./ys + xs.^2 + exp(0.1*ys) - sin(xs);

   [xd, yd] = getBackendVars('diagonal', xv, yv);
   zd = 2*xd - xd.*yd + xd./yd + xd.^2 + exp(0.1*yd) - sin(xd);

   verifyBackendsAgree(t, zs, zd);
end

%--------------------------------------------------------------------------

function testSumAndElementwiseMaxMatchSparseBackend(t)           %#ok<DEFNU>
   % sum() and the elementwise max() branch selection must agree between
   % backends as well, since they take separate code paths (sumJac,
   % plusJac/lMultDiag) from plain arithmetic.
   xv = t.TestData.xv;
   yv = t.TestData.yv;

   [xs, ys] = getBackendVars('sparse', xv, yv);
   zs = sum(max(xs, ys));

   [xd, yd] = getBackendVars('diagonal', xv, yv);
   zd = sum(max(xd, yd));

   verifyBackendsAgree(t, zs, zd);
end

%--------------------------------------------------------------------------

function [a, b] = getBackendVars(kind, av, bv)
   switch kind
      case 'sparse'
         backend = SparseAutoDiffBackend();
      case 'diagonal'
         % Disable MEX acceleration: this suite checks the plain MATLAB
         % diagonal-Jacobian math against the sparse reference backend,
         % and must not depend on a MEX toolchain being available/built.
         backend = DiagonalAutoDiffBackend('useMex', false);
      otherwise
         error('Unknown backend kind %s', kind);
   end
   [a, b] = backend.initVariablesAD(av, bv);
end

%--------------------------------------------------------------------------

function verifyBackendsAgree(t, sparseResult, diagonalResult)
   diagonalResult = diagonalResult.castJacToSparse();

   verifyEqual(t, diagonalResult.val, sparseResult.val, 'AbsTol', 1e-10);

   % NOTE: must copy .jac into a plain cell array before using ':'
   % expansion -- since ADI (and GenericAD) overload subsref, a chained
   % `someADI.jac{:}` only yields the FIRST cs-list element instead of
   % expanding to all of them.
   sjac    = sparseResult.jac;
   djac    = diagonalResult.jac;
   Jsparse = horzcat(sjac{:});
   Jdiag   = horzcat(djac{:});
   verifyEqual(t, full(Jdiag), full(Jsparse), 'AbsTol', 1e-8);
end
