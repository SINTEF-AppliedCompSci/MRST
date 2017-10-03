classdef HandleLinearSolverAD < LinearSolverAD
    % Simple solver for wrapping functions on the form x = fn(A, b);
    properties
        fcn_handle
    end
    
    methods
        function solver = HandleLinearSolverAD(fcn)
            solver.fcn_handle = fcn;
        end
        
       function [result, report] = solveLinearSystem(solver, A, b) %#ok
           report = struct();
           result = solver.fcn_handle(A, b);
       end

    end
end