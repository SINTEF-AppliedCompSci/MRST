classdef CPRSolverAD < LinearSolverAD
    % Solve a problem with a pressure component using constrained a pressure residual method
    %
    % SYNOPSIS:
    %   solver = CPRSolverAD()
    %
    % DESCRIPTION:
    %   Solve a linearized problem with a significant elliptic/pressure
    %   feature via a two stage preconditioner for GMRES. By exposing the
    %   elliptic component as a separate system, a special elliptic solver can
    %   be used to handle the tightly coupled pressure system.
    %
    %   For second stage, ILU(0) is used.
    %
    %
    % OPTIONAL PARAMETERS:
    %   'property' - Set property to value.
    %
    %
    % SEE ALSO:
    %   `BackslashSolverAD`, `LinearSolverAD`, `LinearizedProblem`

    properties
        relativeTolerance % Relative tolerance for elliptic solver
        pressureScaling % Scaling factor applied to pressure equations
        ellipticSolver % LinearSolverAD subclass suitable for the elliptic submatrix.
        diagonalTol % Diagonal tolerance in [0,1].
        ellipticVarName % Name of elliptic-like variable which will be solved using elliptic solver
        trueIMPES % Use true impes decoupling strategy (if supported by model)
        ellipticSign = -1; % Sign to use for pressure system - negative default for compatibility with AGMG
    end
    methods
        function solver = CPRSolverAD(varargin)
            solver = solver@LinearSolverAD();
            
            % Default options
            solver.ellipticSolver = [];
            solver.relativeTolerance = 1e-4;
            solver.pressureScaling = 1/(200*barsa);
            solver.diagonalTol = 1e-2;
            solver.trueIMPES = false;
            solver.ellipticVarName = 'pressure';
            solver.reduceToCell = true;
            
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
        
        function [dx, result, report] = solveLinearProblem(solver, problem0, model)
            % Solve a linearized problem using a constrained pressure
            % residual preconditioner
            timer = tic();
            if ~isempty(problem0.A)
                problem0 = problem0.clearSystem();
                dispif(solver.verbose, 'System already assembled. CPR will re-assemble!');
            end
            % In the case that we have a oil equation, we generally want
            % this to be the first entry as most black oil models are
            % parametrized by oil pressure as the primary variable.
            isOilEq = problem0.indexOfEquationName('oil');
            if any(isOilEq) && ~isOilEq(1)
                indices = 1:numeq(problem0);
                indices(isOilEq) = 1;
                indices(1)   = find(isOilEq);
                problem0 = problem0.reorderEquations(indices);
            end
            
            isPressure = problem0.indexOfPrimaryVariable(solver.ellipticVarName);
            pressureIndex = find(isPressure);
            problem = problem0;
            isCurrent = problem.indexOfType('cell');
            keepNumber0 = solver.keepNumber;
            s = getSampleAD(problem.equations{:});
            % Get pressure indices
            offsets = cumsum([1; s.getNumVars()]);
            pInx = offsets(pressureIndex):offsets(pressureIndex+1)-1;
            if solver.reduceToCell && isempty(solver.keepNumber)
                % Eliminate non-cell variables (well equations etc)
                isCurrent = problem.indexOfType('cell');
                if isa(s, 'GenericAD')
                    % If we are working with block AD, we use the built-in
                    % keepNumber property of the linear solver to perform a
                    % full block Schur complement
                    nk = sum(isCurrent);
                    assert(all(isCurrent(1:nk)) & ~any(isCurrent(nk+1:end)), ...
                        'Cell variables must all combine first in the ordering for this AutodiffBackend.');
                    nv =  s.getNumVars();
                    solver.keepNumber = sum(nv(isCurrent));
                else
                    [problem, eliminated] = solver.reduceToVariable(problem0, isCurrent);
                end
            end
            cellIndex = find(isCurrent);
            
            cellEqNo = numel(cellIndex);
            nCell = problem.getEquationVarNum(1);
            
            % Get and apply scaling
            if solver.trueIMPES
                stype = 'trueIMPES';
            else
                stype = 'simple';
            end
            scale = model.getScalingFactorsCPR(problem, problem.equationNames, stype);
            
            for i = 1:numel(scale)
                if numelValue(scale{i}) > 1 || scale{i} ~= 0
                     problem.equations{i} = problem.equations{i}.*scale{i};
                end
            end
            eqs = problem.equations;
            if isfinite(solver.diagonalTol)
                isElliptic = false(nCell, cellEqNo);
                for i = 1:cellEqNo
                    % Find the derivative of current cell block w.r.t the
                    % pressure
                    pressureJacobi = eqs{i}.jac{isPressure};
                    assert(size(eqs{i}.jac{isPressure}, 2) == nCell);
                    pressureDiag  = diag(pressureJacobi);

                    sod = sum(abs(pressureJacobi), 2) - abs(pressureDiag);
                    % Find "bad" equations, i.e. equations where the current
                    % pressure index does not give a good elliptic jacobian
                    isElliptic(pressureDiag./sod > solver.diagonalTol, i) = true;
                end

                isElliptic(all(isElliptic == 0, 2), pressureIndex) = true;

                bad = ~isElliptic(:, pressureIndex);
                if any(bad)
                    bad = find(bad);
                    % first identify zero diagonal elems of eqs{pressureIndex}
                    isZeroDiag = false(numel(bad), cellEqNo);
                    for i = 1:cellEqNo
                        d = diag(eqs{pressureIndex}.jac{i});
                        isZeroDiag(:,i) = d(bad)==0;
                    end
                    % Switch equations for non-elliptic jacobian components with
                    % some other equation that has an elliptic pressure jacobian
                    [r, c] = find(and(isElliptic(bad, :),~isZeroDiag));
                    sb = numel(bad);
                    if sb == 1
                        % Find gets confused for a single element
                        r = r(:);
                        c = c(:);
                    end
                    replacementInd = accumarray(r, c, [sb, 1], @min);


                    for i = unique(replacementInd)'
                        if i == pressureIndex
                            continue
                        end
                        % Standard switching with some overloaded indexing in
                        % the ad objects
                        isCurrent = bad(replacementInd == i);
                        problem.equations{i}(isCurrent) = eqs{pressureIndex}(isCurrent);
                        problem.equations{pressureIndex}(isCurrent) = eqs{i}(isCurrent);
                    end
                end

                problem.equations{pressureIndex} = isElliptic(:,1).*eqs{pressureIndex};
            else
                isElliptic = true(nCell, cellEqNo);
            end
            for i = 1:cellEqNo
                % Add together all the equations to get a "pressure
                % equation" where we know that all submatrices should be as
                % close as possible to M-matrices (because of switching,
                % which does not actually alter the solution)
                ok = isElliptic(:, i);
                if i == pressureIndex || all(value(scale{i}) == 0)
                    continue
                end
                problem.equations{pressureIndex} = ...
                    problem.equations{pressureIndex} + ok.*eqs{i};
            end
            
            % Solve cell equations
            if solver.pressureScaling ~= 1
                for i = 1:numeq(problem)
                    problem.equations{i}.jac{pressureIndex} =...
                    problem.equations{i}.jac{pressureIndex}./solver.pressureScaling;
                end
            end

            % Set up linear system
            [A, b] = problem.getLinearSystem();
            % Reduce system (if requested)
            [A, b, lsys] = solver.reduceLinearSystem(A, b);
            % ILU0 preconditioner for the non-elliptic part
            [L, U] = ilu(A, struct('type', 'nofill'));
            if isempty(solver.keepNumber)
                Ap = problem.equations{pressureIndex}.jac{pressureIndex};
            else
                Ap = A(pInx, pInx);
            end
            Ap = solver.ellipticSign*Ap;
            % Set up elliptic solver
            solver.ellipticSolver = solver.ellipticSolver.setupSolver(Ap, b(pInx));
            prec = @(r) solver.applyTwoStagePreconditioner(r, A, Ap, L, U, pInx);
            assert(all(isfinite(b)), 'Linear system rhs must have finite entries.');
            t_prep = toc(timer);
            try
                [cprSol, fl, relres, its, resvec] = gmres(A, b, [], solver.tolerance,...
                                                    min(solver.maxIterations, size(A, 1)), prec);
            catch exception
                % Ensure external memory etc is deallocated properly if the
                % solve failed in some spectacular manner
                solver.ellipticSolver = solver.ellipticSolver.cleanupSolver(Ap, b(pInx));
                rethrow(exception)
            end
            t_solve = toc(timer) - t_prep;
            % Clean up elliptic solver
            solver.ellipticSolver = solver.ellipticSolver.cleanupSolver(Ap, b(pInx));

            % Recover eliminated variables on linear level
            cprSol = solver.recoverLinearSystem(cprSol, lsys);
            % Undo pressure scaling
            if solver.pressureScaling ~= 1
                cprSol(pInx) = cprSol(pInx) ./ solver.pressureScaling;
            end
            
            dx = solver.storeIncrements(problem, cprSol);
            
            if solver.reduceToCell && isempty(solver.keepNumber)
                dx = problem.recoverFromSingleVariableType(problem0, dx, eliminated);
            end
            
            solver.keepNumber = keepNumber0;
            t_post = toc(timer) - t_prep - t_solve;
            doPrint = true;
            switch fl
                case 0
                    doPrint = solver.verbose;
                    dispif(doPrint, 'GMRES converged:');
                case 1
                    fprintf('GMRES did not converge. Reached maximum iterations');
                case 2
                    fprintf('GMRES did not converge. Preconditioner was ill-conditioned.');
                case 3
                    fprintf('GMRES stagnated. Unable to reduce residual.');
            end
            dispif(doPrint, ' Final residual: %1.2e after %d iterations (tol: %1.2e) \n', relres, its(2), solver.tolerance);

            if nargout > 1
                result = vertcat(dx{:});
            end
            
            if nargout > 2
                report = solver.getSolveReport(...
                                'Iterations',         its(2), ...
                                'Converged',          fl == 0, ...
                                'SolverTime',         t_solve + t_prep + t_post, ...
                                'LinearSolutionTime', t_solve, ...
                                'PreparationTime',    t_prep, ...
                                'PostProcessTime',    t_post, ...
                                'Residual',           relres);
                report.FlagGMRES = fl;
                if solver.extraReport
                    report.ResidualHistory = resvec;
                end
            end
        end
        
        function [d, sn] = getDescription(solver)
            [tmp, sn_sub] = solver.ellipticSolver.getDescription();
            sn = ['Matlab-CPR-', sn_sub];
            sn = [sn, solver.id];
            d = sprintf(['Matlab implementation of constrained pressure', ...
                        ' residual (CPR) with dynamic row-sum. Elliptic solver: %s.'], sn_sub);
        end
        
        function x = applyTwoStagePreconditioner(solver, r, A, Ap, L, U, pInx)
           es = solver.ellipticSolver;
           
           rp = r(pInx);
           if isfinite(solver.relativeTolerance)
               % Elliptic solver uses absolute tolerance? We scale the
               % problem based on current value.
               n = norm(rp);
               abstol = es.tolerance;
               reltol = solver.relativeTolerance;
               scaling = (abstol/reltol)*(1/n);
           else
               scaling = 1;
           end
           rp = solver.ellipticSign*rp*scaling;
           x = zeros(size(r));
           p = es.solveLinearSystem(Ap, rp)/scaling;
           x(pInx) = p;
           r = r - A*x;
           x = x + U\(L\r);
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
