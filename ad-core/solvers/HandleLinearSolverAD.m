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

        function [d, shortname] = getDescription(solver)
            % Get the description and a short name used for display
            % purposes.
            shortname = 'Handle';
            shortname = [shortname, solver.id];

            d = 'Wrapper for function_handle.';
        end
    end
end

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
