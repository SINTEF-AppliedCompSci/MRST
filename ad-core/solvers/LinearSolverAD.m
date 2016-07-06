classdef LinearSolverAD < handle
%Base class for linear solvers in the AD framework
%
% SYNOPSIS:
%   solver = LinearSolverAD()
%
% DESCRIPTION:
%   Base class for linear solvers. Implement methods for solving linearized
%   problems and adjoints. Also supports setup/cleanup functions
%   before/after use for initialize once-type usage.
%
% REQUIRED PARAMETERS:
%   None
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   See class properties.
%
% NOTE:
%   This class is intended as superclass. It cannot actually solve
%   problems.
%
% SEE ALSO:
%   BackslashSolverAD, CPRSolverAD, LinearizedProblem

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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
   properties
       % Convergence tolerance
       tolerance
       % Max number of iterations used
       maxIterations
       % Enable this to produce additional report output. May use a lot of
       % memory for large problems
       extraReport
       % Verbose output enabler
       verbose
   end
   methods
       function solver = LinearSolverAD(varargin)
           solver.tolerance     = 1e-8;
           solver.maxIterations = 25;
           solver.extraReport   = false;
           solver.verbose       = mrstVerbose();
           solver = merge_options(solver, varargin{:});
           
           assert(solver.maxIterations >= 0);
           assert(solver.tolerance >= 0);
       end
       
       function [result, report] = solveLinearSystem(solver, A, b) %#ok
           % Solve the linear system to a given tolerance
           error('Superclass not meant for direct use')
       end
       
       function [grad, result, report] = solveAdjointProblem(solver, problemPrev,...
                                        problemCurr, adjVec, objective, model) %#ok
           % Solve an adjoint problem.
           problemCurr = problemCurr.assembleSystem();
           
           % Maybe formalize the control variables a bit in the future
           % sometime...
           if iscell(objective)
               objective = objective{:};
           end
           objective = cat(objective);

           rhs = -objective.jac{1}';
           if ~isempty(adjVec)
               problemPrev = problemPrev.assembleSystem();
               rhs = rhs - problemPrev.A'*adjVec;
           end
           A = problemCurr.A';
           [result, report] = solver.solveLinearSystem(A, rhs);
           
           grad = solver.storeIncrements(problemCurr, result);
       end
       
       function [dx, result, report] = solveLinearProblem(solver, problem, model)
           % Solve a linearized problem
           problem = problem.assembleSystem();
           assert(all(isfinite(problem.b)), 'Linear system rhs must have finite entries.');
           
           timer = tic();
           [result, report] = solver.solveLinearSystem(problem.A, problem.b); 
           report.SolverTime = toc(timer);
           
           dx = solver.storeIncrements(problem, result);
       end
       
       function dx = storeIncrements(solver, problem, result) %#ok
           % Extract the results from a vector into a cell array with one
           % entry per primary variable in the linearized problem.
           
           % Find first index corresponding to ADI equation
           ix = find(cellfun(@(x) isa(x, 'ADI'), problem.equations), 1);
           % Calculate positions in newton increment
           numVars = cellfun(@(x) size(x, 2), problem.equations{ix}.jac)';
           cumVars = cumsum(numVars);
           ii = [[1;cumVars(1:end-1)+1], cumVars];
           
           eqn = size(ii,1);
           dx = cell(eqn,1);
           for i = 1:eqn
               dx{i} = result(ii(i,1):ii(i,2));
           end
       end
       
       function solver = setupSolver(solver, A, b, varargin) %#ok 
           % Run setup on a solver for a given system
           
           % Dummy function run before a set of linear problems with
           % different right hand sides
       end
       
       function  solver = cleanupSolver(solver, A, b, varargin) %#ok 
           % Clean up solver after use (if needed)
           
           % Dummy function run after a set of linear problems with
           % different right hand sides 
           
           % For overloading when for example calling multigrid solvers as
           % a preconditioner
       end
       
        function [dx, result, report] = solveCellReducedLinearProblem(solver, problem, model)
            % Reduce a problem to cell-variables, solve and then recover
            % the eliminated variables afterwards.
            
            % Eliminate non-cell variables (well equations etc)
            keep = problem.indexOfType('cell');
            
            [problem, eliminated] = solver.reduceToVariable(problem, keep);
            
            % Solve a linearized problem
            problem = problem.assembleSystem();
            
            timer = tic();
            [result, report] = solver.solveLinearSystem(problem.A, problem.b);
            report.SolverTime = toc(timer);
            
            dxCell = solver.storeIncrements(problem, result);
            
            dx = solver.recoverResult(dxCell, eliminated, keep);
        end
        
        function [problem, eliminated] = reduceToVariable(solver, problem, keep)
            remove = find(~keep);
            
            problem = problem.clearSystem();
            
            eliminated = cell(numel(remove), 1);
            elimNames = problem.equationNames(remove);
            
            for i = 1:numel(remove)
                [problem, eliminated{i}] = problem.eliminateVariable(elimNames{i});
            end
        end
        
        function dx = recoverResult(solver, dxElim, eliminatedEqs, keep)
            kept = find(keep);
            left = find(~keep);
            keptEqNo = numel(kept);
            
            % Find number of variables
            nP = numel(keep);
            
            % Set up storage for all variables, including those we
            % eliminated previously
            dx = cell(nP, 1);
            
            % Recover non-cell variables
            recovered = false(nP, 1);
            recovered(kept) = true;
            
            % Put the recovered variables into place
            dx(recovered) = dxElim;
            
            for i = numel(eliminatedEqs):-1:1
                pos = left(i);
                dVal = recoverVars(eliminatedEqs{i}, keptEqNo + 1, dx(recovered));
                dx{pos} = dVal;
                recovered(pos) = true;
            end
        end
   end
end
