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
       doApplyScalingCPR % true / false
       useSYMRCMOrdering % true / false
       pressureScaling % scaling factor for pressure - automatically determined
       diagonalTol = 0.2; % tolerance if strategy ends with _drs
       couplingTol = 0.02; % tolerance for drs
       decoupling = 'trueIMPES'; % trueimpes, quasiimpes, none
       strategy = 'mrst'; % mrst, mrst_drs, amgcl, amgcl_drs
   end
   methods
       function solver = AMGCL_CPRSolverAD(varargin)
            require linearsolvers
            solver = solver@AMGCLSolverAD();
            solver.doApplyScalingCPR = true;
            solver.reduceToCell = true;
            solver.tolerance    = 1e-6;
            solver.useSYMRCMOrdering = false;

            [solver, extra] = merge_options(solver, varargin{:});
            solver.amgcl_setup = getAMGCLMexStruct(extra{:});
            solver.amgcl_setup.solver_id = 2;
       end

       function [result, report] = solveLinearSystem(solver, A, b, varargin)
           [result, report] = solver.callAMGCL_MEX(A, b, solver.amgcl_setup.solver_id, varargin{:});
       end

       function [dx, result, report] = solveLinearProblem(solver, problem, model)
           problem = solver.prepareProblemCPR(problem, model);
           [dx, result, report] = solveLinearProblem@LinearSolverAD(solver, problem, model);
       end

       function [dx, result, report] = solveAdjointProblem(solver, problemPrev,problemCurr, adjVec, objective, model)
           if ~isempty(problemPrev)
                problemPrev = solver.prepareProblemCPR(problemPrev, model);
           end
           problemCurr = solver.prepareProblemCPR(problemCurr, model);
           [dx, result, report] = solveAdjointProblem@LinearSolverAD(solver, problemPrev,problemCurr, adjVec, objective, model);
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
       
       function setSRelaxation(solver, varargin)
           solver.setParameterGroup('relaxation', 's_relaxation', varargin{:});
       end

       function problem = prepareProblemCPR(solver, problem, model)
           n = model.G.cells.num;
           solver.pressureScaling = mean(value(problem.state.pressure));
           setup = solver.amgcl_setup;
           if setup.update_sprecond || setup.update_ptransfer
               % Reuse AMG hierarchy
               solver.reuseMode = 2;
               if problem.iterationNo == 1
                   resetAMGCL();
               end
           else
               solver.reuseMode = 1;
           end
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
                if ~isempty(problem.A)
                    problem = problem.clearSystem();
                    dispif(solver.verbose, 'System already assembled. CPR will re-assemble!');
                end
               scale = model.getScalingFactorsCPR(problem, problem.equationNames, solver.decoupling);
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

           if isempty(solver.keepNumber)
               if solver.reduceToCell
                   % Will be reduced to ncell by block_size system
                   ndof = n*m;
               else
                   % We have no idea and should check
                   problem = problem.assembleSystem();
                   ndof = size(problem.A, 1);
                   if solver.amgcl_setup.active_rows == 0
                       % Only the first n*m entries are cell-wise
                       % variables, tell the solver this
                       solver.amgcl_setup.active_rows = n*m;
                   end
               end
           else
               ndof = solver.keepNumber;
           end

           if isempty(solver.variableOrdering) || numel(solver.variableOrdering) ~= ndof
               if solver.useSYMRCMOrdering
                   sym_ordering = getGridSYMRCMOrdering(model.G);
               else
                   sym_ordering = [];
               end
               ordering = getCellMajorReordering(n, m, 'ndof', ndof, 'cell_ordering', sym_ordering);
               solver.variableOrdering = ordering;
               if isempty(solver.equationOrdering) || numel(solver.equationOrdering) ~= ndof
                   solver.equationOrdering = ordering;
               end
           end

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

        function [A, b, scaling, x0] = applyScaling(solver, A, b, x0)
            if nargin == 3
                x0 = [];
            end
            bz = solver.amgcl_setup.block_size;
            nc = solver.amgcl_setup.cell_size;
            psub = (1:bz:(nc*bz - bz + 1))';

            if solver.applyLeftDiagonalScaling || solver.applyRightDiagonalScaling
                [A, b, scaling, x0] = applyScaling@LinearSolverAD(solver, A, b, x0);
            elseif numel(solver.pressureScaling) ~= 1 || solver.pressureScaling ~= 1
                n = size(A, 1);
                d = ones(n, 1);
                d(~isfinite(d)) = 1;
                d(psub) = solver.pressureScaling;
                I = (1:n)';
                M = sparse(I, I, d, n, n);
                A = A*M;
                if ~isempty(x0)
                    x0(psub) = x0(psub)./d(psub);
                end
                scaling.M = M;
            end

            useMRST = strcmpi(solver.strategy, 'mrst');
            useMRST_drs = strcmpi(solver.strategy, 'mrst_drs');

            if useMRST || useMRST_drs
                w = getScalingInternalCPR(solver, A, b);

                ncv = bz*nc;
                ndof = size(b, 1);

                I = rldecode((1:bz:ncv)', bz);
                J = (1:ncv)';

                if 0
                    D = sparse(I, J, w, ndof, ndof);
                    tmp = D*A;
                    btmp = D*b;
                    A(psub, :) = tmp(psub, :);
                    b(psub) = btmp(psub);
                else
                    Id = J;
                    Id(1:bz:(ncv-bz+1)) = [];
                    Id = [Id; ((ncv+1):ndof)'];
                    D = sparse([I; Id], ...
                               [J; Id], ...
                               [w(1:ncv); ones(numel(Id), 1)], ndof, ndof);
                    A = D*A;
                    b = D*b;
                end
                w_override = zeros(ncv, 1);
                w_override(psub) = 1;
                solver.amgcl_setup.drs_row_weights = w_override;
            elseif strcmpi(solver.strategy, 'amgcl_drs')
                solver.amgcl_setup.drs_row_weights = getScalingInternalCPR(solver, A, b);
            end
        end

        function x = undoScaling(solver, x, scaling) %#ok
            % Undo effects of scaling applied to linear system
            if isfield(scaling, 'M') && ~isempty(scaling.M)
                x = scaling.M*x;
            end
        end

        function M = getDiagonalInverse(solver, A)
            % Reciprocal of diagonal matrix
            if solver.applyRightDiagonalScaling
                sz = size(A);
                assert(sz(1) == sz(2), 'Matrix must be square!');
                bz = solver.amgcl_setup.block_size;
                nc = solver.amgcl_setup.cell_size;
                n = sz(1);
                d = 1./abs(diag(A));
                d(~isfinite(d)) = 1;
                psub = 1:bz:(nc*bz - bz + 1);
                d(psub) = solver.pressureScaling;
                I = (1:n)';
                M = sparse(I, I, d, n, n);
            else
                M = getDiagonalInverse@LinearSolverAD(solver, A);
            end
        end
   end
end

function [w, ndof] = getScalingInternalCPR(solver, A, b)
    bz = solver.amgcl_setup.block_size;
    nc = solver.amgcl_setup.cell_size;
    ndof = bz*nc;
    p_inx = 1:bz:(ndof - bz + 1);
    switch lower(solver.decoupling)
        case 'quasiimpes'
            [ii, jj, vv] = find(A);
            blockNoI = floor((ii-1)/bz)+1;
            blockNoJ = floor((jj-1)/bz)+1;
            keep = blockNoJ >= blockNoI & blockNoJ < (blockNoI+1)  & ii <= ndof & jj <= ndof;
            if 1
                % Weights are determined by constant in rhs
                Ap = sparse(jj(keep), ii(keep), vv(keep), ndof, ndof);
            else
                % Weights sum up to unity
                isp = false(ndof, 1);
                isp(p_inx) = true;
                I = jj(keep);
                J = ii(keep);
                V = vv(keep);

                V(isp(I) & ~isp(J)) = 1;
                V(isp(I) & isp(J)) = 1;

                Ap = sparse(I, J, V, ndof, ndof);
            end
            q = zeros(ndof, 1);
            q(p_inx) = 1;
            w = Ap\q;

            w = [w; ones(size(A, 1)-ndof, 1)];

        otherwise
            w = ones(size(b, 1), 1);
    end

    if strcmpi(solver.strategy, 'mrst_drs')
        if not(strcmpi(solver.decoupling, 'quasiimpes'))
            [ii, jj, vv] = find(A);
        end
        % Dynamic row sum strategy by Gries et al, SPE-163608-PA.
        isp = false(numel(ii), 1);
        isp(p_inx) = true;
        blockNo = ceil(ii./bz);
        blockConn = ceil(jj./bz);
        isdp = isp(jj);
        isdiag = blockConn == blockNo & isdp;

        is_offdiag_p = isdp & ~isdiag;

        pd = zeros(ndof, 1);
        pd(ii(isdiag)) = vv(isdiag);
        cellno = ceil((1:ndof)'/bz);
        % Check diagonal dominance
        if isfinite(solver.diagonalTol)
            e_dd = solver.diagonalTol;
            sum_offdiag = accumarray(ii(is_offdiag_p), abs(vv(is_offdiag_p)), [ndof, 1]);

            ok_dd = pd >= sum_offdiag*e_dd;
        else
            ok_dd = true(ndof, 1);
        end
        % Check for very weak coupling
        if isfinite(solver.couplingTol)
            e_ps = solver.couplingTol;
            is_other_p = isdp & blockConn ~= blockNo;
            sum_other_blocks = accumarray(ii(is_other_p), abs(vv(is_other_p)), [ndof, 1]);
            ok_ps = sum_other_blocks >= pd*e_ps;
        else
            ok_ps = true(ndof, 1);
        end

        ok = ok_dd & ok_ps;
        % Check if all cells are bad
        ok_count = accumarray(cellno, ok);
        bad = ok_count == 0;
        if any(bad)
            % All equations are somehow not diagonally
            % dominant, pick the first entry.
            ok((find(bad)-1)*bz+1) = true;
        end
        % Always set first entry to true for the time being. We should
        % probably switch equations if these are not ok, but this is not
        % implemented.
        ok(p_inx) = true;

        w(~ok) = 0;
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
