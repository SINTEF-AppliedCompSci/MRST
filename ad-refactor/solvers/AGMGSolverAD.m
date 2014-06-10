classdef AGMGSolverAD < linearSolverAD
    % Base class for a nonlinear solver
   properties
       
   end
   methods
       function solver = AGMGSolverAD(varargin)
            solver = solver@linearSolverAD(varargin{:});
       end
       
       function result = solveLinearSystem(solver, A, b)
           % Solve the linear system to a given tolerance
           result = agmg(A, b, [], solver.tolerance, solver.maxiterations, [], [], 2);
       end
       
       function solver = setupSolver(solver, A, b, varargin) %#ok 
           % Run setup on a solver for a given system
           agmg(A,[],[],[],[],[],[], 1);
       end
       
       function  solver = cleanupSolver(solver, A, b, varargin) %#ok 
           % Clean up solver after use (if needed)
           agmg(A,[],[],[],[],[],[], -1);
       end
   end
end