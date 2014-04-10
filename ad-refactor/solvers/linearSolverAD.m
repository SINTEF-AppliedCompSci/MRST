classdef linearSolverAD
    % Base class for a nonlinear solver
   properties
       tolerance
       maxiterations
   end
   methods
       function solver = linearSolverAD(varargin)
           opt = struct('tolerance',     1e-8, ...
                        'maxiterations', 25);
            
           solver.tolerance = opt.tolerance;
           solver.maxiterations = opt.maxiterations;
       end
       
       function result = solveLinearSystem(solver, A, b) %#ok
           error('Superclass not meant for direct use')
       end
       
       function [dx, result] = solveLinearProblem(solver, problem)
           problem = problem.assembleSystem();
           result = solver.solveLinearSystem(problem.A, problem.b); 
           
           % Store the increments seperately - copy/paste from SolveEqsADI
           numVars = cellfun(@numval, problem.equations)';
           cumVars = cumsum(numVars);
           ii = [[1;cumVars(1:end-1)+1], cumVars];
           
           eqn = size(ii,1);
           dx = cell(eqn,1);
           for i = 1:eqn
               dx{i} = result(ii(i,1):ii(i,2));
           end
       end
       
       function setupSolver(solver, A, b, varargin) %#ok 
           % Dummy function run before a set of linear problems with
           % different right hand sides
       end
       
       function cleanupSolver(solver, A, b, varargin) %#ok 
           % Dummy function run after a set of linear problems with
           % different right hand sides 
           
           % For overloading when for example calling multigrid solvers as
           % a preconditioner
       end
   end
end