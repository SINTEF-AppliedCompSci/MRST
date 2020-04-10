classdef AMGCL_CPRSolverBlockAD < AMGCLSolverBlockAD
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
       doApplyScalingCPR % true / false
       pressureScaling % scaling factor for pressure - automatically determined
       diagonalTol = 0.2; % tolerance if strategy ends with _drs
       couplingTol = 0.02; % tolerance for drs
       decoupling = 'trueIMPES'; % trueimpes, quasiimpes, none
       strategy = 'mrst'; % mrst, mrst_drs, amgcl, amgcl_drs
   end
   methods
       function solver = AMGCL_CPRSolverBlockAD(varargin)
            require linearsolvers
            solver = solver@AMGCLSolverBlockAD();
            solver.doApplyScalingCPR = false;
            solver.reduceToCell = true;
            solver.tolerance    = 1e-6;

            [solver, extra] = merge_options(solver, varargin{:});
            solver.amgcl_setup = getAMGCLMexStruct(extra{:});
            % solver.amgcl_setup.solver_id = 2;
            solver.amgcl_id = 2;
       end

       function [result, report] = solveLinearSystem(solver, A, b)
           [result, report] = solver.callAMGCL_MEX(A, b, solver.amgcl_setup.solver_id);
       end

       function [dx, result, report] = solveLinearProblem(solver, problem, model)
           problem = solver.prepareProblemCPR(problem, model);
           [dx, result, report] = solveLinearProblem@AMGCLSolverBlockAD(solver, problem, model);
       end
       
       function setSRelaxation(solver, varargin)
           solver.setParameterGroup('relaxation', 's_relaxation', varargin{:});
       end

       function problem = prepareProblemCPR(solver, problem, model)
           n = model.G.cells.num;
           solver.pressureScaling = mean(value(problem.state.pressure));

           if solver.amgcl_setup.block_size == 0
               % Solver has not been told about block size, try to compute
               % it from what we are given.
               s = getSampleAD(problem.equations{:});
               nv = s.getNumVars();
               isCell = problem.indexOfType('cell');
               solver.amgcl_setup.block_size = sum(nv(isCell)/n);
           end
           solver.amgcl_setup.cell_num = model.G.cells.num;
           solver.amgcl_setup.cell_size = n;


           % Get and apply scaling
           if solver.doApplyScalingCPR && strcmpi(solver.decoupling, 'trueimpes')
               assert(false)
                if ~isempty(problem.A)
                    problem = problem.clearSystem();
                    dispif(solver.verbose, 'System already assembled. CPR will re-assemble!');
                end
               scale = model.getScalingFactorsCPR(problem, problem.equationNames, solver);
               % Solver will take the sum for us, we just weight each
               % equation. Note: This is not the entirely correct way
               % of doing this, as solver could do this by itself.
               for i = 1:numel(scale)
                   if ~strcmpi(problem.types{i}, 'cell')
                       continue
                   end
                   ds = value(scale{i});
                   if (numel(ds) > 1 || any(ds ~= 0))
                       problem.equations{i} = problem.equations{i}.*scale{i};
                   end
               end
           end
           m = solver.amgcl_setup.block_size;
           assert(m > 0);

           switch lower(solver.strategy)
               case {'mrst', 'mrst_drs'}
                   solver.amgcl_setup.use_drs = true;
                   solver.amgcl_setup.drs_eps_ps = -1e8;
                   solver.amgcl_setup.drs_eps_dd = -1e8;
               case 'amgcl'
                   solver.amgcl_setup.use_drs = false;
               case 'amgcl_drs'
                   solver.amgcl_setup.use_drs = true;
                   if isfinite(solver.couplingTol)
                       solver.amgcl_setup.drs_eps_ps = solver.couplingTol;
                   else
                       solver.amgcl_setup.drs_eps_ps = -1e8;
                   end
                   if isfinite(solver.diagonalTol)
                       solver.amgcl_setup.drs_eps_dd = solver.diagonalTol;
                   else
                       solver.amgcl_setup.drs_eps_dd = -1e8;
                   end
               otherwise
                   error('Unknown CPR strategy %s', solver.strategy);
           end
       end
       
       function [d, sn] = getDescription(solver)
           sn = 'AMGCL-CPR';
           if solver.amgcl_setup.cpr_blocksolver
               sn = [sn, '-block'];
           end
           sn = [sn, solver.id];
           prm = {'solver', 'preconditioner', 'relaxation'};
           if solver.amgcl_setup.preconditioner == 1
               prm{end+1} = 'coarsening';
           end
           tmp = cell(1, numel(prm)+1);
           for i = 1:numel(prm)
               s = prm{i};
               tmp{i} = solver.getFeatureDescription(s);
           end
           tmp{end} = solver.getFeatureDescription('relaxation', 's_relaxation');
           d = [sprintf('AMGCL constrained-pressure-residual (CPR) solver. Configuration:\n'), ...
               sprintf('\t%s\n', tmp{:})];
       end
   end
end


%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
