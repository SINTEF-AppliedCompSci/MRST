classdef mldivideSolverAD < linearSolverAD
   methods
       function result = solveLinearSystem(solver, A, b) %#ok
           result = A\b;
       end
   end
end