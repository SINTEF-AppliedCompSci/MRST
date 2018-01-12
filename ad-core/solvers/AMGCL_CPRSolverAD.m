classdef AMGCL_CPRSolverAD < AMGCLSolverAD
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
       doApplyScalingCPR
       trueIMPES % Use true impes decoupling strategy (if supported by model)
       
   end
   methods
       function solver = AMGCL_CPRSolverAD(varargin)
            require linearsolvers
            solver = solver@AMGCLSolverAD();
            solver.trueIMPES    = true;
            solver.doApplyScalingCPR = true;
            solver.reduceToCell = true;
            solver.tolerance    = 1e-6;
            
            [solver, extra] = merge_options(solver, varargin{:});
            solver.amgcl_setup = getAMGCLMexStruct(extra{:});
       end
       
       function [result, report] = solveLinearSystem(solver, A, b)
            report = struct();
            [result, error_estimate] = amgcl_matlab_cpr(A', b, solver.amgcl_setup);
            if error_estimate > solver.tolerance
                warning('Solver did not converge to specified tolerance of %g. Reported residual estimate was %g', solver.tolerance, error_estimate);
            end
       end
       
        function [dx, result, report] = solveLinearProblem(solver, problem, model)
            % Get and apply scaling
            if solver.doApplyScalingCPR
                scale = model.getScalingFactorsCPR(problem, problem.equationNames, solver);
                e = 0;
                for i = 1:numel(scale)
                    if ~strcmpi(problem.types{i}, 'cell')
                        continue
                    end
                    if (numel(scale{i}) > 1 || scale{i} ~= 0)
                         e = e + problem.equations{i}.*scale{i};
                    end
                end
                problem.equations{1} = e;
            end
            
            n = model.G.cells.num;
            m = solver.amgcl_setup.block_size;
            assert(m > 0);
            ndof = n*m;
            if isempty(solver.variableOrdering) || numel(solver.variableOrdering) ~= ndof
                ordering = getCellMajorReordering(n, m, ndof);
                solver.variableOrdering = ordering;
                if isempty(solver.equationOrdering) || numel(solver.equationOrdering) ~= ndof
                    solver.equationOrdering = ordering;
                end
            end
            
            [dx, result, report] = solveLinearProblem@LinearSolverAD(solver, problem, model);
        end
        
       function setSRelaxation(solver, v)
           solver.amgcl_setup.s_relaxation = translateOptionsAMGCL('s_relaxation', v);
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
