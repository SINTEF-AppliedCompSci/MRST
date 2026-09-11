function tests = LinearizedProblemTest
%Test suite for LinearizedProblem (Newton-loop linear system assembly)
%
% DESCRIPTION:
%   `LinearizedProblem` assembles the linear system solved at every
%   Newton iteration and performs the block-Gaussian elimination used
%   to remove well/control equations before calling a linear solver,
%   but had no dedicated unit test -- only indirect exercise via
%   full-simulation regression tests. This suite checks
%   `assembleSystem` against a hand-assembled Jacobian/residual for a
%   small two-equation ADI system, a handful of query helpers, and the
%   Schur-complement elimination performed by `eliminateVariable`
%   against a hand-solved reduced system.
%
% SEE ALSO:
%   LinearizedProblem, NonLinearSolver, functiontests

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
end

%--------------------------------------------------------------------------

function testAssembleSystemMatchesHandAssembledJacobianAndResidual(t) %#ok<DEFNU>
   % For two residual equations eq1 = 2*x + 3*y - 10 and
   % eq2 = x - y + 1 over a 2-cell vector x, y, the assembled linear
   % system must be the block matrix [[2*I,3*I];[I,-I]] and the
   % right-hand side must be -[eq1.val; eq2.val].
   x0 = [1; 2]; y0 = [3; 4];
   [X, Y] = initVariablesADI(x0, y0);
   eq1 = 2*X + 3*Y - 10;
   eq2 = X - Y + 1;

   problem = LinearizedProblem({eq1, eq2}, {'cell', 'cell'}, {'eq1', 'eq2'}, ...
                                {'x', 'y'}, struct());
   problem = problem.assembleSystem();

   I2 = eye(2);
   expectedA = [2*I2, 3*I2; I2, -I2];
   expectedb = -[eq1.val; eq2.val];

   verifyEqual(t, full(problem.A), expectedA, 'AbsTol', 1e-10);
   verifyEqual(t, problem.b, expectedb, 'AbsTol', 1e-10);
end

%--------------------------------------------------------------------------

function testNormAndQueryHelpers(t)                              %#ok<DEFNU>
   % norm() must give the per-equation vector norm of each residual;
   % numeq/countOfType/indexOfPrimaryVariable/indexOfEquationName must
   % report consistently with how the problem was constructed.
   x0 = [1; 2]; y0 = [3; 4];
   [X, Y] = initVariablesADI(x0, y0);
   eq1 = 2*X + 3*Y - 10;
   eq2 = X - Y + 1;

   problem = LinearizedProblem({eq1, eq2}, {'cell', 'cell'}, {'eq1', 'eq2'}, ...
                                {'x', 'y'}, struct());

   verifyEqual(t, problem.norm(), [norm(eq1.val), norm(eq2.val)], 'AbsTol', 1e-10);
   verifyEqual(t, numeq(problem), 2);
   verifyEqual(t, problem.countOfType('cell'), 2);
   verifyEqual(t, problem.indexOfPrimaryVariable('y'), [false, true]);
   verifyEqual(t, problem.indexOfEquationName('eq2'), [false, true]);
   verifyEqual(t, problem.getEquationVarNum(), [2; 2]);
end

%--------------------------------------------------------------------------

function testEliminateVariablePerformsCorrectSchurComplement(t)  %#ok<DEFNU>
   % Eliminating equation/variable 'X' from
   %   eqX = X + 2*Y - 10   (dX=1, dY=2)
   %   eqY = 3*X - Y - 4    (dX=3, dY=-1)
   % via the block-Gaussian formula
   %   eqY.jac{Y} -= eqY.jac{X} * (eqX.jac{X} \ eqX.jac{Y})
   %   eqY.val    -= eqY.jac{X} * (eqX.jac{X} \ eqX.val)
   % must give a reduced single-equation problem with
   %   dEqY/dY = -1 - 3*(1\2) = -7,  val = 4 - 3*(1\9) = -23
   % (using x0=5, y0=7, so eqX.val=9, eqY.val=4).
   x0 = 5; y0 = 7;
   [X, Y] = initVariablesADI(x0, y0);
   eqX = X + 2*Y - 10;
   eqY = 3*X - Y - 4;

   problem = LinearizedProblem({eqX, eqY}, {'well', 'cell'}, {'X', 'Y'}, ...
                                {'X', 'Y'}, struct());
   [reduced, eliminated] = problem.eliminateVariable('X');

   verifyEqual(t, eliminated.val, 9);
   verifyEqual(t, reduced.primaryVariables, {'Y'});
   verifyEqual(t, reduced.equationNames, {'Y'});
   verifyEqual(t, reduced.types, {'cell'});
   verifyEqual(t, full(reduced.equations{1}.jac{1}), -7, 'AbsTol', 1e-10);
   verifyEqual(t, reduced.equations{1}.val, -23, 'AbsTol', 1e-10);
end
