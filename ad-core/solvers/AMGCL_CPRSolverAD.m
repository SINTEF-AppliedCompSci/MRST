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
       useSYMRCMOrdering
       pressureScaling
   end
   methods
       function solver = AMGCL_CPRSolverAD(varargin)
            require linearsolvers
            solver = solver@AMGCLSolverAD();
            solver.trueIMPES    = false;
            solver.doApplyScalingCPR = true;
            solver.reduceToCell = true;
            solver.tolerance    = 1e-6;
            solver.useSYMRCMOrdering = true;
            
            [solver, extra] = merge_options(solver, varargin{:});
            solver.amgcl_setup = getAMGCLMexStruct(extra{:});
       end
       
       function [result, report] = solveLinearSystem(solver, A, b)
           [result, report] = solver.callAMGCL_MEX(A, b, 2);
%            bz = solver.amgcl_setup.block_size;
%            nc = solver.amgcl_setup.cell_size;
%            
%            result(1:bz:(nc*bz)) = result(1:bz:(nc*bz))/(1000*barsa);
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

       function setSRelaxation(solver, v)
           solver.amgcl_setup.s_relaxation = translateOptionsAMGCL('relaxation', v);
       end
       
       function problem = prepareProblemCPR(solver, problem, model)
           n = model.G.cells.num;
           solver.pressureScaling = mean(problem.state.pressure);
           if solver.amgcl_setup.block_size == 0
               % Solver has not been told about block size, try to compute
               % it from what we are given.
               s = getSampleAD(problem.equations{:});
               nv = s.getNumVars();
               isCell = problem.indexOfType('cell');
               solver.amgcl_setup.block_size = sum(nv(isCell)/n);
           end
           solver.amgcl_setup.cell_num = model.G.cells.num;
                for ph = 1:size(problem.state.s, 2)
%                     problem.equations{ph} = problem.equations{ph}.*(problem.dt./model.operators.pv);
                end
            if isa(model, 'ThreePhaseBlackOilModel') && false
                sw = model.getProp(problem.state, 'sw');
                bad = double(sw) <= 1e-8;
                if any(bad)
                    sum(bad)
                    for ph = 2:size(problem.state.s, 2)
%                         problem.equations{1}(bad) = problem.equations{1}(bad) + problem.equations{ph}(bad);
                    end
                end
           end
           solver.amgcl_setup.cell_size = n;

           
           % Get and apply scaling
           if solver.doApplyScalingCPR  && false
               if 1
                   scale = model.getScalingFactorsCPR(problem, problem.equationNames, solver);
                   % Solver will take the sum for us, we just weight each
                   % equation. Note: This is not the entirely correct way
                   % of doing this, as solver could do this by itself.
                   for i = 1:numel(scale)
                       if ~strcmpi(problem.types{i}, 'cell')
                           continue
                       end
                       ds = double(scale{i});
                       if (numel(ds) > 1 || any(ds ~= 0))
                           problem.equations{i} = problem.equations{i}.*scale{i};
                       end
                   end
               else
               end
           end
           m = solver.amgcl_setup.block_size;
           assert(m > 0);
           
           if isempty(solver.keepNumber)
               if solver.reduceToCell
                   % Will be reduced to ncell by block_size syste,
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
       end
       
        function [A, b, scaling] = applyScaling(solver, A, b)
            bz = solver.amgcl_setup.block_size;
            nc = solver.amgcl_setup.cell_size;
            psub = (1:bz:(nc*bz - bz + 1))';

            if solver.applyLeftDiagonalScaling || solver.applyRightDiagonalScaling
                [A, b, scaling] = applyScaling@LinearSolverAD(solver, A, b);
            elseif solver.pressureScaling ~= 1
                n = size(A, 1);
                d = ones(n, 1);
                d(~isfinite(d)) = 1;
                d(psub) = solver.pressureScaling;
                I = (1:n)';
                M = sparse(I, I, d, n, n);
                A = A*M;
                scaling.M = M;
            end
            
            if solver.amgcl_setup.use_drs
                [w, ncv] = getScalingInternalCPR(solver, A, b);
                
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
                solver.amgcl_setup.drs_eps_dd = -1e8;
                solver.amgcl_setup.drs_eps_dd = -1e8;
                
                w_override = zeros(ncv, 1);
                w_override(psub) = 1;
                solver.amgcl_setup.drs_row_weights = w_override;
            end
            if 0
                w = getScalingInternalCPR(solver, A, b);
                ix = (1:numel(w))';
                D = sparse(ix, ix, w);
                A = D*A;
                b = D*b;
            end
        end
        
        function x = undoScaling(solver, x, scaling)
            x = undoScaling@LinearSolverAD(solver, x, scaling);
        end
        
        function M = getDiagonalInverse(solver, A)
            % Reciprocal of diagonal matrix
            if solver.applyRightDiagonalScaling
                sz = size(A);
                assert(sz(1) == sz(2), 'Matrix must be square!');
                bz = solver.amgcl_setup.block_size;
                nc = solver.amgcl_setup.cell_size;
                n = sz(1);
%                 d = 1./abs(diag(A));
                d = 1./diag(A);
%                 d = ones(n, 1);
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
    [ii, jj, vv] = find(A);
    
    bz = solver.amgcl_setup.block_size;
    nc = solver.amgcl_setup.cell_size;
    ndof = bz*nc;
    
    p_inx = 1:bz:(ndof - bz + 1);

    
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
    
    
%     if solver.applyRightDiagonalScaling
        q(p_inx) = 1;
%     else
%         q(p_inx) = 1e-8;
%     end
    w = Ap\q;
%     
%     tmp = reshape(w, 2, [])';
%     tmp = 1e-3*tmp./sum(tmp, 2);
%     tmp = tmp';
%     w = tmp(:);
    
    w = [w; ones(size(A, 1)-ndof, 1)];
end

%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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