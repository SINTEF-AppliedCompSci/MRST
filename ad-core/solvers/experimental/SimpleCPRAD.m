classdef SimpleCPRAD < CPRSolverAD
% Solve a problem with a pressure component using constrained a pressure residual method
%
% SYNOPSIS:
%   solver = CPRSolverAD()
%
% DESCRIPTION:
%   Solve a linearized problem with a significant elliptic/pressure
%   component via a two stage preconditioner for GMRES. By exposing the
%   elliptic component as a seperate system, a special elliptic solver can
%   be used to handle the highly connected components.
%
%   For second stage, ILU(0) is used.
%
% REQUIRED PARAMETERS:
%   None
%
% OPTIONAL PARAMETERS:
%   See class properties.
%
%
% SEE ALSO:
%   BackslashSolverAD, LinearSolverAD, LinearizedProblem

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
    properties
        eliminateAll
    end
    methods
        function solver = SimpleCPRAD(varargin)
            solver = solver@CPRSolverAD();
            solver.ellipticSolver = [];
            solver.eliminateAll = true;
            % Default options
            solver = merge_options(solver, varargin{:});
                        
            if isempty(solver.ellipticSolver)
                solver.ellipticSolver = BackslashSolverAD();
            else
                assert(isa(solver.ellipticSolver, 'LinearSolverAD'));
            end
        end
        
        function [result, report] = solveLinearSystem(solver, A, b) %#ok
            error('Not supported directly - this is a preconditioner')
        end
        
        function [dx, result, report] = solveLinearProblem(solver, problem, model)
            % Solve a linearized problem using a constrained pressure
            % residual preconditioner
            timer = tic();
            
            % Get and apply scaling
            scale = model.getScalingFactorsCPR(problem, problem.equationNames, solver);
            
            for i = 1:numel(scale)
                if numel(scale{i}) > 1 || scale{i} ~= 0
                     problem.equations{i} = problem.equations{i}.*scale{i};
                end
            end
            nPh = model.water + model.oil + model.gas;
            for i = 2:nPh
                problem.equations{1} = problem.equations{1} + problem.equations{i};
            end

            % Set up linear system
            [A, b] = problem.getLinearSystem();
            
            n = size(A, 1);
            nc = model.G.cells.num;
            

            if solver.eliminateAll
                keepNum = (model.water + model.gas + model.oil)*nc;
                [reorderEq, reorderVar] = deal([]);
            else
                eqCounts = cellfun(@(x) numel(double(x)), problem.equations)';
                reorderEq = nan(n, 1);
                neq = numeq(problem);
                move = false(neq, 1);
                for i = 1:neq
                    move(i) = any(strcmpi(problem.equationNames{i}, ...
                        {'segClosure', 'pDropSeg', 'waterwells', 'oilwells', 'gaswells', 'closurewells'}));
                end
                isMove = rldecode(move, eqCounts);
                nmove = sum(eqCounts.*move);
                reorderEq(isMove) = (n-nmove+1:n)';%(1:nmove)';
                reorderEq(~isMove) = (1:n-nmove)';
                reorderVar = reorderEq;
                keepNum = n - nmove;
            end
%             neq = numeq(problem); 
%             eqs = problem.equations; 
%             for ix = 1:neq, 
%                 if full(min(abs(diag(eqs{ix}.jac{ix})))) == 0
%                     fprintf('Pair %s, %s has zero diagonal\n', problem.primaryVariables{ix},...
%                         problem.equationNames{ix});
%                 end
%             end
            
            
            %tic()
            p = eliminate(A, b, keepNum, reorderEq, reorderVar);
            %toc();
            
            A = p.A;
            b = p.b;
            
            % ILU0 preconditioner for the non-elliptic part
            [L, U] = ilu(A, struct('type', 'nofill'));
            
            pInx = 1;
            pInx = (1:nc).' + nc*(pInx-1);
            Ap = A(pInx, pInx);
            % We have only cell variables present, and these will have
            % offsets of cellnum long each
            
            % Set up elliptic solver
            solver.ellipticSolver = solver.ellipticSolver.setupSolver(Ap, b(pInx));
            ellipSolve = @(b) solver.ellipticSolver.solveLinearSystem(Ap, b);

            prec = @(r) applyTwoStagePreconditioner(r, A, L, U, pInx, ellipSolve);
            assert(all(isfinite(b)), 'Linear system rhs must have finite entries.');
            try
                [result, fl, relres, its, resvec] = gmres(A, b, [], solver.relativeTolerance,...
                                                    min(solver.maxIterations, size(A, 1)), prec);
            catch exception
                % Ensure external memory etc is deallocated properly if the
                % solve failed in some spectacular manner
                solver.ellipticSolver = solver.ellipticSolver.cleanupSolver(Ap, b(pInx));
                rethrow(exception)
            end
            % Clean up elliptic solver
            solver.ellipticSolver = solver.ellipticSolver.cleanupSolver(Ap, b(pInx));
            
            result = recover(result, p);
            
            dx = solver.storeIncrements(problem, result);
            
            % Recover stuff
            solvetime = toc(timer);
            
            if solver.verbose
                switch fl
                    case 0
                        fprintf('GMRES converged:');
                    case 1
                        fprintf('GMRES did not converge. Reached maximum iterations');
                    case 2
                        fprintf('GMRES did not converge. Preconditioner was ill-conditioned.');
                    case 3
                        fprintf('GMRES stagnated. Unable to reduce residual.');
                end
                fprintf(' Final residual: %1.2e after %d iterations (tol: %1.2e) \n', relres, its(2), solver.relativeTolerance);
            end
            
            if nargout > 1
                result = vertcat(dx{:});
            end
            
            if nargout > 2
                report = struct('IterationsGMRES', its(2), ...
                                'FlagGMRES',       fl, ...
                                'SolverTime',      solvetime, ...
                                'FinalResidual',   relres);
                if solver.extraReport
                    report.ResidualHistory = resvec;
                end
            end
        end
        
    end
end

function x = applyTwoStagePreconditioner(r, A, L, U, pInx, ellipticSolver)
   x = zeros(size(r));
   x(pInx) = ellipticSolver(r(pInx));

   r = r - A*x;
   x = x + U\(L\r);
end
function p = eliminate(A, b, keepNum, reorderEq, reorderVar)
   n = size(A, 1);
   
   start = 1:keepNum;
   
   if ~isempty(reorderVar)
       A = A(reorderEq, reorderVar);
   end
   [ix, jx, vx] = find(A);
   
   p = struct();
   if 1
       keep = false(n, 1);
       keep(start) = true;
       nk = keepNum;
       
       keepRow = keep(ix);
       keepCol = keep(jx);
       kb = keepRow & keepCol;
       p.B = sparse(ix(kb), jx(kb), vx(kb), nk, nk);
       
       kc = keepRow & ~keepCol;
       p.C = sparse(ix(kc), jx(kc) - nk, vx(kc), nk, n - nk);
       
       kd = ~keepRow & keepCol;
       p.D = sparse(ix(kd) - nk, jx(kd), vx(kd), n - nk, nk);
       
       ke = ~keepRow & ~keepCol;
       p.E = sparse(ix(ke) - nk, jx(ke) - nk, vx(ke), n - nk, n - nk);
       p.f = b(keep);
       p.h = b(~keep);
   else
       p.B = A(start, start);

       p.C = A(start, stop);
       p.D = A(stop, start);
       p.E = A(stop, stop);

       p.f = b(start);
       p.h = b(stop);
   end
   [Le, Ue] = lu(p.E);
   p.A = p.B - p.C*(Ue\(Le\p.D));
   p.b = p.f - p.C*(Ue\(Le\p.h));

%    p.A = p.B - p.C*(p.E\p.D);
%    p.b = p.f - p.C*(p.E\p.h);
   
   p.reorder = reorderVar;
end

function result = recover(result, p)
   s = p.E\(p.h - p.D*result);
   result = [result; s];
   if ~isempty(p.reorder)
        backmap = zeros(numel(result), 1);
        backmap(p.reorder) = (1:numel(result))';
        result = result(backmap);
   end
end
