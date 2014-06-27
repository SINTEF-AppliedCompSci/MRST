classdef NoOpSolverAD < linearSolverAD
    % Base class for a nonlinear solver
   properties
       
   end
   methods
       
       function [result, report] = solveLinearSystem(solver, A, b)
           result = zeros(size(b));
           report = struct('Converged', false);
       end
   end
end