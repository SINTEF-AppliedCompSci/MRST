classdef mldivideSolverAD < linearSolverAD
   methods
       function [result, report] = solveLinearSystem(solver, A, b) %#ok
           result = A\b;
           % Nothing to report
           report = struct();
       end
   end
end