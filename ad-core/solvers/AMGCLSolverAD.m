classdef AMGCLSolverAD < LinearSolverAD
    % Linear solver that calls external compiled multigrid solver
    %
    % SYNOPSIS:
    %   solver = AMGCLSolverAD()
    %
    % DESCRIPTION:
    %    AD-interface for the AMGCL interface.
    %
    % NOTE:
    %    This solver requires AMGCL to be installed and working.
    %
    % SEE ALSO:
    %   `BackslashSolverAD`

   properties
        amgcl_setup
   end
   methods
       function solver = AMGCLSolverAD(varargin)
            require linearsolvers
            solver = solver@LinearSolverAD();
            solver.reduceToCell = true;
            solver.tolerance = 1e-6;
            
            [solver, extra] = merge_options(solver, varargin{:});
            solver.amgcl_setup = getAMGCLMexStruct(extra{:});

       end
       
       function [result, report] = solveLinearSystem(solver, A, b)
            report = struct();
            [result, error_estimate] = amgcl_matlab(A', b, solver.amgcl_setup);
            if error_estimate > solver.tolerance
                warning('Solver did not converge to specified tolerance of %g. Reported residual estimate was %g', solver.tolerance, error_estimate);
            end
       end
       
       function setCoarsening(solver, v)
           solver.amgcl_setup.coarsening = translateOptionsAMGCL('coarsening', v);
       end
       
       function setRelaxation(solver, v)
           solver.amgcl_setup.relaxation = translateOptionsAMGCL('relaxation', v);
       end
       
       function setPreconditioner(solver, v)
           solver.amgcl_setup.preconditioner = translateOptionsAMGCL('preconditioner', v);
       end
   end
end

%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

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
