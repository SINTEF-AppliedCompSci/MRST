function tests = ADITest
%Test suite for the ADI automatic-differentiation class (ADI, initVariablesADI)
%
% SEE ALSO:
%   ADI, initVariablesADI, functiontests

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

function testConstructorGivesIdentityJacobians(t)                %#ok<DEFNU>
   % initVariablesADI must give each variable an identity Jacobian with
   % respect to itself, and a correctly-sized all-zero Jacobian with
   % respect to every other variable in the same call.
   x0 = [1; 2; 3];
   y0 = [4; 5];

   [X, Y] = initVariablesADI(x0, y0);

   verifyEqual(t, X.val, x0);
   verifyEqual(t, Y.val, y0);
   verifyEqual(t, X.jac{1}, speye(3));
   verifyEqual(t, full(X.jac{2}), zeros(3, 2));
   verifyEqual(t, full(Y.jac{1}), zeros(2, 3));
   verifyEqual(t, Y.jac{2}, speye(2));
end

%--------------------------------------------------------------------------

function testAdditionAndSubtraction(t)                           %#ok<DEFNU>
   % d(x+y)/dx = I, d(x+y)/dy = I; d(x-y)/dx = I, d(x-y)/dy = -I.
   x0 = [1; 2; 3];
   y0 = [4; 5; 6];
   [X, Y] = initVariablesADI(x0, y0);

   S = X + Y;
   verifyEqual(t, S.val, x0 + y0);
   verifyEqual(t, S.jac{1}, speye(3));
   verifyEqual(t, S.jac{2}, speye(3));

   D = X - Y;
   verifyEqual(t, D.val, x0 - y0);
   verifyEqual(t, D.jac{1}, speye(3));
   verifyEqual(t, full(D.jac{2}), full(-speye(3)));
end

%--------------------------------------------------------------------------

function testElementwiseMultiplicationAndDivision(t)             %#ok<DEFNU>
   % d(x.*y)/dx = diag(y), d(x.*y)/dy = diag(x).
   % d(x./y)/dx = diag(1./y), d(x./y)/dy = diag(-x./y.^2).
   x0 = [1; 2; 3];
   y0 = [4; 5; 6];
   [X, Y] = initVariablesADI(x0, y0);

   P = X .* Y;
   verifyEqual(t, P.val, x0 .* y0);
   verifyEqual(t, full(P.jac{1}), diag(y0));
   verifyEqual(t, full(P.jac{2}), diag(x0));

   Q = X ./ Y;
   verifyEqual(t, Q.val, x0 ./ y0, 'AbsTol', 1e-12);
   verifyEqual(t, full(Q.jac{1}), diag(1./y0), 'AbsTol', 1e-12);
   verifyEqual(t, full(Q.jac{2}), diag(-x0./y0.^2), 'AbsTol', 1e-12);
end

%--------------------------------------------------------------------------

function testPower(t)                                            %#ok<DEFNU>
   % d(x.^2)/dx = diag(2*x), and the power must carry a correctly-sized
   % all-zero Jacobian with respect to unrelated variables.
   x0 = [1; 2; 3];
   y0 = [4; 5; 6];
   [X, Y] = initVariablesADI(x0, y0); %#ok<ASGLU>

   Z = X.^2;
   verifyEqual(t, Z.val, x0.^2);
   verifyEqual(t, full(Z.jac{1}), diag(2*x0));
   verifyEqual(t, full(Z.jac{2}), zeros(3));
end

%--------------------------------------------------------------------------

function testElementaryFunctions(t)                               %#ok<DEFNU>
   % Derivatives of exp, log, sin and cos must match their closed forms.
   x0 = [0.5; 1; 1.5];
   X  = initVariablesADI(x0);

   E = exp(X);
   verifyEqual(t, E.val, exp(x0));
   verifyEqual(t, full(E.jac{1}), diag(exp(x0)), 'AbsTol', 1e-12);

   L = log(X);
   verifyEqual(t, L.val, log(x0));
   verifyEqual(t, full(L.jac{1}), diag(1./x0), 'AbsTol', 1e-12);

   S = sin(X);
   verifyEqual(t, S.val, sin(x0));
   verifyEqual(t, full(S.jac{1}), diag(cos(x0)), 'AbsTol', 1e-12);

   C = cos(X);
   verifyEqual(t, C.val, cos(x0));
   verifyEqual(t, full(C.jac{1}), diag(-sin(x0)), 'AbsTol', 1e-12);
end

%--------------------------------------------------------------------------

function testIndexingSubsref(t)                                  %#ok<DEFNU>
   % Indexing an ADI object must slice both value and Jacobian rows
   % consistently.
   x0 = [1; 2; 3; 4];
   X  = initVariablesADI(x0);

   sub = X([2, 4]);
   I = speye(4);
   verifyEqual(t, sub.val, x0([2, 4]));
   verifyEqual(t, sub.jac{1}, I([2, 4], :));
end

%--------------------------------------------------------------------------

function testIndexingSubsasgn(t)                                 %#ok<DEFNU>
   % Assigning a plain (non-AD) value into a subset of an ADI object must
   % overwrite that subset's value and zero out the corresponding
   % Jacobian rows.
   x0 = [1; 2; 3; 4];
   X  = initVariablesADI(x0);

   X([1, 3]) = [10; 30];
   I = speye(4);
   verifyEqual(t, X.val, [10; 2; 30; 4]);
   verifyEqual(t, full(X.jac{1}([1, 3], :)), zeros(2, 4));
   verifyEqual(t, full(X.jac{1}([2, 4], :)), full(I([2, 4], :)));
end

%--------------------------------------------------------------------------

function testVertcat(t)                                          %#ok<DEFNU>
   % Vertical concatenation must stack both values and Jacobian rows.
   x0 = [1; 2];
   X  = initVariablesADI(x0);

   Z = [X; X];
   verifyEqual(t, Z.val, [x0; x0]);
   verifyEqual(t, full(Z.jac{1}), [full(speye(2)); full(speye(2))]);
end

%--------------------------------------------------------------------------

function testSumAndSingleArgMax(t)                                %#ok<DEFNU>
   % sum(X).jac must be the row-sum of X.jac; max(X) (single argument)
   % must pick out the value and Jacobian row of the maximizing element.
   x0 = [3; 1; 2];
   X  = initVariablesADI(x0);

   S = sum(X);
   verifyEqual(t, S.val, sum(x0));
   verifyEqual(t, full(S.jac{1}), ones(1, 3));

   M = max(X);
   verifyEqual(t, M.val, 3);
   verifyEqual(t, full(M.jac{1}), [1, 0, 0]);
end

%--------------------------------------------------------------------------

function testElementwiseMaxAndMin(t)                             %#ok<DEFNU>
   % Elementwise max/min must select, per element, the value and
   % Jacobian row of whichever operand is larger/smaller. At a tie
   % (index 3, where x0 == y0), both max and min consistently pick the
   % second argument's branch (max compares with a strict '>', and min
   % is implemented as -max(-u,-v), which preserves that same tie-break).
   x0 = [1; 5; 3];
   y0 = [4; 2; 3];
   [X, Y] = initVariablesADI(x0, y0);

   Mx = max(X, Y);
   verifyEqual(t, Mx.val, [4; 5; 3]);
   verifyEqual(t, full(Mx.jac{1}), diag([0, 1, 0]));
   verifyEqual(t, full(Mx.jac{2}), diag([1, 0, 1]));

   Mn = min(X, Y);
   verifyEqual(t, Mn.val, [1; 2; 3]);
   verifyEqual(t, full(Mn.jac{1}), diag([1, 0, 0]));
   verifyEqual(t, full(Mn.jac{2}), diag([0, 1, 1]));
end

%--------------------------------------------------------------------------

function testCompositeExpressionMatchesFiniteDifferenceJacobian(t) %#ok<DEFNU>
   % A composite, purely elementwise expression must have an ADI
   % Jacobian diagonal matching a central finite-difference
   % approximation. Since every output element depends only on the
   % input element of the same index, this is a strong end-to-end check
   % of the chain rule through +, .*, ./, .^, sin and exp together.
   f = @(x) sin(x).*exp(-x) + x.^3./(1 + x);

   x0 = linspace(0.3, 2.1, 6)';
   X  = initVariablesADI(x0);
   F  = f(X);

   h  = 1e-6;
   fd = zeros(size(x0));
   for i = 1:numel(x0)
      xp    = x0; xp(i) = xp(i) + h;
      xm    = x0; xm(i) = xm(i) - h;
      fdVal = (f(xp) - f(xm)) ./ (2*h);
      fd(i) = fdVal(i);
   end

   verifyEqual(t, F.val, f(x0), 'AbsTol', 1e-12);
   verifyEqual(t, diag(full(F.jac{1})), fd, 'AbsTol', 1e-6);
end
