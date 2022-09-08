classdef LinearSolverAD < handle
    % Base class for linear solvers in the AD framework
    %
    % SYNOPSIS:
    %   solver = LinearSolverAD()
    %
    % DESCRIPTION:
    %   Base class for linear solvers. Implement methods for solving linearized
    %   problems and adjoints. Also supports setup/cleanup functions
    %   before/after use for initialize once-type usage.
    %
    % REQUIRED PARAMETERS:
    %   None
    %
    % OPTIONAL PARAMETERS:
    %   'property' - Set property to value.
    %
    % NOTE:
    %   This class is intended as superclass. It cannot actually solve
    %   problems.
    %
    % SEE ALSO:
    %   `BackslashSolverAD`, `CPRSolverAD`, `LinearizedProblem`
    
    properties
        tolerance % Linear solver tolerance
        maxIterations % Max number of iterations used
        extraReport % Enable this to produce additional report output
        % May use a lot of memory for large problems
        verbose % Verbose output enabler
        replaceNaN % Boolean indicating if the solver should replace NaN in the results
        replaceInf % Boolean indicating if the solver should replace Inf in the results
        replacementNaN % If replaceNaN is enabled, this is the value that will be inserted
        replacementInf % If replaceInf is enabled, this is the value that will be inserted
        reduceToCell % Reduce to per-cell system before solving
        applyLeftDiagonalScaling % Apply left diagonal scaling before solving
        applyRightDiagonalScaling % Apply right diagonal scaling before solving
        keepNumber % If set, linear solver will reduce the system to the first keepNumber entries
        useSparseReduction % If true, sparse indexing will be used with keepNumber option
        variableOrdering % Variable ordering to be used for linear solver. Row vector of equal length to the size of the linear system.
        equationOrdering % Equation ordering to be used for linear solver. Row vector of equal length to the size of the linear system.
        id = ''; % Short text string identifying the specific solver. Appended to the short name (see getDescription)
        initialGuessFun = [];
    end

    methods
        function solver = LinearSolverAD(varargin)
            solver.tolerance       = 1e-8;
            solver.maxIterations   = 25;
            solver.extraReport     = false;
            solver.verbose         = mrstVerbose();
            solver.replaceNaN      = false;
            solver.replaceInf      = false;
            solver.replacementNaN  = 0;
            solver.replacementInf  = 0;
            solver.reduceToCell    = false;
            solver.useSparseReduction = false;
            solver.applyLeftDiagonalScaling = false;
            solver.applyRightDiagonalScaling = false;
            solver.variableOrdering = [];
            solver.equationOrdering = [];
            
            solver = merge_options(solver, varargin{:});
            
            assert(solver.maxIterations >= 0);
            assert(solver.tolerance >= 0);
        end
        
        function [result, report] = solveLinearSystem(solver, A, b, varargin) %#ok
            % Solve the linear system to a given tolerance
            report = solver.getSolveReport();
            error('Superclass not meant for direct use')
        end
        
        function report = getSolveReport(solver, varargin) %#ok
            report = struct('Iterations', 0, ... % Number of iterations (if iterative)
                            'Residual',   0, ... % Final residual
                            'SolverTime', 0, ... % Total time in solver
                            'LinearSolutionTime', 0, ... % Time spent solving system
                            'PreparationTime', 0, ... % Schur complement, scaling, ...
                            'PostProcessTime', 0, ... % Recovery, undo scaling, ...
                            'Converged', true); % Bool indicating convergence
            report = merge_options_relaxed(report, varargin);
        end
        
        function [grad, result, report] = solveAdjointProblem(solver, problemPrev,...
                problemCurr, adjVec, objective, model, varargin) %#ok
            
            opt = struct('scalePressure',  false, ...
                         'colIx',            nan,   ...
                         'equationScaling',   []);
            % For the adjoint problem, the scaling for the rows of the matrix to be inverted
            % corresponds to pressure and, say, saturation. This bad scaling
            % triggers typically a warning from Matlab which says that the
            % matrix ill-conditioned. To prevent that, we can rescale the
            % pressure. Note that, even if the condition number of the matrix
            % is unchanged by transposition, no warning is typically
            % triggered for the forward problem.
            % The implementation is not necessarily robust with variable/equation reordering
            % or equation reduction, that is why it is now only implemented as a switch.
            opt = merge_options(opt, varargin{:});
            
            % Solve an adjoint problem.
            timer = tic();
            problemCurr = problemCurr.assembleSystem();
            
            % Maybe formalize the control variables a bit in the future
            % sometime...
            if iscell(objective)
                objective = objective{:};
            end
            objective = combineEquations(objective);
            if isempty(objective) || isempty(objective.val)
                % no explicit dependence on this time-step
                b = [];
            else
                assert(isa(objective, 'ADI'), 'Objective function was not of type ADI.');
                b = -objective.jac{1}';
            end

            if ~isempty(adjVec)
                % handle pre-allocated zero adjVec for last step:
                if isempty(problemPrev)
                    problemPrev = problemCurr;
                end
                % hack
                ix = find(cellfun(@(x)isa(x, 'GenericAD'), problemPrev.equations));
                if any(ix)
                    if ~isempty(b)
                        mismatch = numel(b) - sum(problemPrev.equations{ix(1)}.numVars);
                        if mismatch ~= 0
                            % adjust size of jacobians
                            nd = numel(problemPrev.equations) - numel(ix);
                            assert(mod(mismatch, nd)==0, 'Unable to resolve jacobian mismatch');
                            for k = 1:numel(ix)
                                problemPrev.equations{ix(k)}.jac{2}.dim(1) = ...
                                    problemPrev.equations{ix(k)}.jac{2}.dim(1) + mismatch/nd;
                            end
                        end
                    end
                    nad = cellfun(@numelValue, problemPrev.equations(ix));
                    ix = find(cellfun(@(x)isa(x, 'double'), problemPrev.equations));
                    nd  = cellfun(@numel, problemPrev.equations(ix));
                    mismatch = numel(adjVec) - sum(nad) - sum(nd);
                    if mismatch ~= 0
                        assert(mod(mismatch, numel(ix))==0, 'Unable to resolve jacobian mismatch');
                        for k  =1:numel(ix)
                            problemPrev.equations{ix(k)} = zeros(nd(k) + mismatch/numel(ix),1);
                        end
                    end
                end
                        
                problemPrev = problemPrev.assembleSystem();
                if isempty(b)
                    b = - problemPrev.A'*adjVec;
                else
                    if ~isfinite(opt.colIx)
                        % standard scalar objective
                        b = b - problemPrev.A'*adjVec;
                    else
                        % Matrix right-hand-side 
                        if max(opt.colIx) > size(adjVec,2)
                            warning('Insufficient size of Lagrange muliplier matrix');
                        end
                        tmp = b;
                        b = - problemPrev.A'*adjVec;
                        b(:, opt.colIx) = b(:, opt.colIx) + tmp;
                    end
                end
            end
            A = problemCurr.A;
            if isempty(b)
                % should only happen in for scalar objectives
                b = zeros(size(A,1),1);
            end
            b = full(b);
            if opt.scalePressure
                stateCurr = problemCurr.state;
                p = model.getProps(stateCurr, 'pressure');
                pmax = max(p);
                np = numel(p);
                [nz, nz] = size(A);
                D = speye(nz, nz);
                ind = sub2ind(size(A), 1 : np, 1 : np);
                D(ind) = pmax*ones(np, 1);
                t_prepare = toc(timer);
                % Solve system
                [result, report] = solver.solveLinearSystem(D'*A, D'*b);
                t_solve = toc(timer) - t_prepare;
            else
                % Reduce system (if requested)
                [A, b, lsys] = solver.reduceLinearSystemAdjoint(A, b);
                % Reorder linear system
                [A, b, ordering] = solver.reorderLinearSystem(A, b);
                % Apply scaling
                [A, b, scaling] = solver.applyScalingAdjoint(A, b);
                % Apply transpose
                A = A';
                t_prepare = toc(timer);
                % Solve system
                [result, report] = solver.solveLinearSystem(A, b);
                t_solve = toc(timer) - t_prepare;
                % Undo scaling
                result = solver.undoScalingAdjoint(result, scaling);
                % Permute system back
                result = solver.deorderLinearSystemAdjoint(result, ordering);
                % Recover eliminated variables on linear level
                result = solver.recoverLinearSystemAdjoint(result, lsys);
            end
            eqScale = opt.equationScaling;
            if ~isempty(eqScale)
                if iscell(eqScale)
                    eqScale = vertcat(eqScale{:});
                end
                result = eqScale.*result;
            end
            report.SolverTime = toc(timer);
            report.LinearSolutionTime = t_solve;
            report.PreparationTime = t_prepare;
            report.PostProcessTime = report.SolverTime - t_solve - t_prepare;
            grad = solver.storeIncrements(problemCurr, result);
        end
        
        function [dx, result, report] = solveLinearProblem(solver, problem, model)
            % Solve a linearized problem
            timer = tic();
            lsys = [];
            eliminated = {};
            keepNumber0 = solver.keepNumber;
            initialGuess = solver.getInitialGuess(problem);
            if solver.reduceToCell && isempty(solver.keepNumber)
                % Eliminate non-cell variables (well equations etc)
                s = getSampleAD(problem.equations{:});
                keep = problem.indexOfType('cell');
                if ~all(keep)
                    if isa(s, 'GenericAD')
                        % If we are working with block AD, we use the built-in
                        % keepNumber property of the linear solver to perform a
                        % full block Schur complement
                        nk = sum(keep);
                        assert(all(keep(1:nk)) & ~any(keep(nk+1:end)), ...
                            'Cell variables must all combine first in the ordering for this AutodiffBackend.');
                        if 1
                            % In-place Schur complement
                            ngroups = numel(s.offsets)-1;
                            varno = rldecode((1:ngroups)', diff(s.offsets));

                            keepEq = problem.equations(keep);
                            elimEq = problem.equations(~keep);

                            keepVar = unique(varno(keep));
                            elimVar = setdiff(varno, keepVar);

                            B_eq = keepEq;
                            C_eq = B_eq;
                            for i = 1:numel(B_eq)
                                B_eq{i}.jac = B_eq{i}.jac(keepVar);
                                C_eq{i}.jac = C_eq{i}.jac(elimVar);
                            end
                            D_eq = elimEq;
                            E_eq = D_eq;
                            for i = 1:numel(D_eq)
                                D_eq{i}.jac = D_eq{i}.jac(keepVar);
                                E_eq{i}.jac = E_eq{i}.jac(elimVar);
                            end
                            B_eq = combineEquations(B_eq{:});
                            C_eq = combineEquations(C_eq{:});
                            D_eq = combineEquations(D_eq{:});
                            E_eq = combineEquations(E_eq{:});
                            lsys = struct('B', B_eq.jac{1}, ...
                                          'C', C_eq.jac{1}, ...
                                          'D', D_eq.jac{1}, ...
                                          'E', E_eq.jac{1},...
                                          'f', -B_eq.val, ...
                                          'h', -D_eq.val, ...
                                          'E_L', [], ...
                                          'E_U', []);
                            [lsys.E_L, lsys.E_U] = lu(lsys.E);
                            problem.A = lsys.B - lsys.C*(lsys.E_U\(lsys.E_L\lsys.D));
                            problem.b = lsys.f - lsys.C*(lsys.E_U\(lsys.E_L\lsys.h));
                        else
                            nv =  s.getNumVars();
                            solver.keepNumber = sum(nv(keep));
                        end
                    else
                        [problem, eliminated] = solver.reduceToVariable(problem, keep);
                    end
                    if ~isempty(initialGuess{1})
                        initialGuess = initialGuess(keep);
                    end
                end
            end
            problem = problem.assembleSystem();
            assert(all(isfinite(problem.b)), 'Linear system rhs must have finite entries.');

            % Get linearized system
            [A, b] = problem.getLinearSystem();
            x0     = vertcat(initialGuess{:});
            % Reduce system (if not already done)
            if isempty(lsys)
                [A, b, lsys, x0] = solver.reduceLinearSystem(A, b, false, x0);
            end
            % Reorder linear system
            [A, b, ordering, x0] = solver.reorderLinearSystem(A, b, [], x0);
            % Apply scaling
            [A, b, scaling, x0] = solver.applyScaling(A, b, x0);

            t_prepare = toc(timer);
            % Solve the system
            [result, report] = solver.solveLinearSystem(A, b, x0);
            t_solve = toc(timer) - t_prepare;
            % Undo scaling
            result = solver.undoScaling(result, scaling);
            % Permute system back
            result = solver.deorderLinearSystem(result, ordering);
            % Recover eliminated variables on linear level
            result = solver.recoverLinearSystem(result, lsys);
            
            [result, report] = problem.processResultAfterSolve(result, report);
            report.SolverTime = toc(timer);
            report.LinearSolutionTime = t_solve;
            report.PreparationTime = t_prepare;
            report.PostProcessTime = report.SolverTime - t_solve - t_prepare;
            if solver.replaceNaN
                result(isnan(result)) = solver.replacementNaN;
            end
            if solver.replaceInf
                result(isinf(result)) = solver.replacementInf;
            end
            dx = solver.storeIncrements(problem, result);
            if ~isempty(eliminated)
                dx = solver.recoverResult(dx, eliminated, keep);
            end
            solver.keepNumber = keepNumber0;
        end
        
        function dx = storeIncrements(solver, problem, result) %#ok
            % Extract the results from a vector into a cell array with one
            % entry per primary variable in the linearized problem.
            
            % Find first index corresponding to ADI equation
            ix = find(cellfun(@(x) isa(x, 'ADI'), problem.equations), 1);
            if isempty(ix)
                % No ADI equations, we return empty increments
                assert(isempty(result), ...
                    'No ADI equations returned. Unable to map increments to variables.');
                dx = {};
                return
            end
            % Calculate positions in newton increment
            numVars = problem.equations{ix}.getNumVars();
            cumVars = cumsum(numVars);
            ii = [[1;cumVars(1:end-1)+1], cumVars];
            
            eqn = size(ii,1);
            dx = cell(eqn,1);
            for i = 1:eqn
                dx{i} = result(ii(i,1):ii(i,2), :);
            end
        end
        
        function [A, b, sys, x0] = reduceLinearSystem(solver, A, b, isAdjoint, x0)
            % Perform Schur complement reduction of linear system
            if nargin < 4
                isAdjoint = false;
            end
            if nargin < 5
                x0 = [];
            end
            if solver.useSparseReduction
                method = 'sparse';
            else
                method = 'subset';
            end
            sys = splitMatrixForReduction(A, b, solver.keepNumber, method, true);
            if ~isempty(sys.B)
                A = sys.B - sys.C*(sys.E_U\(sys.E_L\sys.D));
                if isAdjoint
                    % We are solving the transpose system
                    b = sys.f - (sys.D')*((sys.E')\sys.h);
                else
                    % We are solving the untransposed system
                    b = sys.f - sys.C*(sys.E_U\(sys.E_L\sys.h));
                end
                if ~isempty(x0)
                    x0 = x0(1:solver.keepNumber);
                end
            end
        end
        
        function [A, b, sys] = reduceLinearSystemAdjoint(solver, A, b)
            [A, b, sys] = reduceLinearSystem(solver, A, b, true);
        end

        function x = recoverLinearSystem(solver, x, sys)
            % Recover eliminated variables
            if ~isempty(sys.E_U)
                s = sys.E_U\(sys.E_L\(sys.h - sys.D*x));
                x = [x; s];
            end
        end
        
        function x = recoverLinearSystemAdjoint(solver, x, sys)
            % Recover eliminated variables
            if ~isempty(sys.E)
                s = (sys.E')\(sys.h - sys.C'*x);
                x = [x; s];
            end
        end

        function [A, b, scaling, x0] = applyScaling(solver, A, b, x0)
            % Apply left or right diagonal scaling
            if nargin < 4
                x0 = [];
            end
            scaling = struct();
            applyLeft = solver.applyLeftDiagonalScaling;
            applyRight = solver.applyRightDiagonalScaling;
            if ~applyLeft && ~applyRight
                return
            end
            M = solver.getDiagonalInverse(A);
            if solver.applyLeftDiagonalScaling
                assert(~applyRight, 'Cannot both apply left and right diagonal scaling');
                A = M*A;
                b = M*b;
            else
                A  = A*M;
                if ~isempty(x0)
                    x0 = M\x0;
                end
            end
            scaling.M = M;
        end
        
        function [A, b, scaling, x0] = applyScalingAdjoint(solver, A, b, x0)
            % Apply left or right diagonal scaling
            if nargin < 4
                x0 = [];
            end
            scaling = struct();
            applyLeft = solver.applyLeftDiagonalScaling;
            applyRight = solver.applyRightDiagonalScaling;
            if ~applyLeft && ~applyRight
                return
            end
            M = solver.getDiagonalInverse(A);
            if solver.applyLeftDiagonalScaling
                assert(~applyRight, 'Cannot both apply left and right diagonal scaling');
                A = M*A;
                if ~isempty(x0)
                    x0 = M'\x0;
                end
            else
                A = A*M;
                b = M'*b;
            end
            scaling.M = M;
        end
        
        function x = undoScaling(solver, x, scaling)
            % Undo effects of scaling applied to linear system
            if solver.applyRightDiagonalScaling
                x = scaling.M*x;
            end
        end
        
        function x = undoScalingAdjoint(solver, x, scaling)
            % Undo effects of scaling applied to linear system (adjoint
            % version)
            if solver.applyLeftDiagonalScaling
                x = scaling.M'*x;
            end
        end

        function x = preconditionerInverse(solver, M, x)
            % Apply a preconditioner. Either a handle or a matrix.
            if isempty(M)
                return
            end
            
            if isa(M, 'function_handle')
                % We got a function handle for the inverse
                x = M(x);
            else
                % We got a matrix
                x = M\x;
            end
        end
        
        function [A, b, order, x0] = reorderLinearSystem(solver, A, b, order, x0)
            if nargin < 5
                x0 = [];
            end
            vo = solver.variableOrdering;
            eo = solver.equationOrdering;
            hasVar = ~isempty(vo);
            hasEq = ~isempty(eo);
            n = size(A, 1);
            if hasVar
                if isa(vo, 'function_handle')
                    vo = vo(A, b);
                end
                nv = numel(vo);
                if nv < n
                    vo = [vo; (nv+1:n)'];
                elseif nv > n
                    vo = vo(1:n);
                end
            end
            if hasEq
                if isa(eo, 'function_handle')
                    eo = eo(A, b);
                elseif isnan(eo)
                    eo = vo;
                end
                ne = numel(eo);
                if ne < n
                    eo = [eo; (ne+1:n)'];
                elseif ne > n
                    eo = eo(1:n);
                end
            end
            if hasVar && hasEq
                A  = A(eo, vo);
                b  = b(eo);
            elseif hasVar
                A  = A(:, vo);
            elseif hasEq
                A = A(eo, :);
                b = b(eo);
            end
            if ~isempty(x0) && hasVar
                x0 = x0(vo);
            end
            order = struct('variableOrdering', vo, 'equationOrdering', eo);
        end
        
        function x = deorderLinearSystem(solver, x, order)
            if nargin < 3
                tmp = solver;
            else
                tmp = order;
            end
            if ~isempty(solver.variableOrdering)
                x(tmp.variableOrdering) = x(1:numel(tmp.variableOrdering));
            end
        end

        function x = deorderLinearSystemAdjoint(solver, x, order)
            if nargin < 3
                tmp = solver;
            else
                tmp = order;
            end
            if ~isempty(solver.equationOrdering)
                x(tmp.equationOrdering) = x(1:numel(tmp.equationOrdering));
            end
        end

        function solver = setupSolver(solver, A, b, varargin) %#ok
            % Run setup on a solver for a given system
            
            % Dummy function run before a set of linear problems with
            % different right hand sides
        end
        
        function  solver = cleanupSolver(solver, A, b, varargin) %#ok
            % Clean up solver after use (if needed)
            
            % Dummy function run after a set of linear problems with
            % different right hand sides
            
            % For overloading when for example calling multigrid solvers as
            % a preconditioner
        end
        
        function [dx, result, report] = solveCellReducedLinearProblem(solver, problem, model)
            % Reduce a problem to cell-variables, solve and then recover
            % the eliminated variables afterwards.
            
            % Eliminate non-cell variables (well equations etc)
            keep = problem.indexOfType('cell');
            
            [problem, eliminated] = solver.reduceToVariable(problem, keep);
            
            % Solve a linearized problem
            problem = problem.assembleSystem();
            
            timer = tic();
            [result, report] = solver.solveLinearSystem(problem.A, problem.b);
            report.SolverTime = toc(timer);
            
            dxCell = solver.storeIncrements(problem, result);
            
            dx = solver.recoverResult(dxCell, eliminated, keep);
        end
        
        function [problem, eliminated] = reduceToVariable(solver, problem, keep)
            remove = find(~keep);
            
            problem = problem.clearSystem();
            
            eliminated = cell(numel(remove), 1);
            elimNames = problem.equationNames(remove);
            
            for i = 1:numel(remove)
                [problem, eliminated{i}] = problem.eliminateVariable(elimNames{i});
            end
        end
        
        function dx = recoverResult(solver, dxElim, eliminatedEqs, keep)
            kept = find(keep);
            left = find(~keep);
            keptEqNo = numel(kept);
            
            % Find number of variables
            nP = numel(keep);
            
            % Set up storage for all variables, including those we
            % eliminated previously
            dx = cell(nP, 1);
            
            % Recover non-cell variables
            recovered = false(nP, 1);
            recovered(kept) = true;
            
            % Put the recovered variables into place
            dx(recovered) = dxElim;
            
            for i = numel(eliminatedEqs):-1:1
                pos = left(i);
                dVal = recoverVars(eliminatedEqs{i}, keptEqNo + 1, dx(recovered));
                dx{pos} = dVal;
                recovered(pos) = true;
            end
        end
        
        function M = getDiagonalInverse(solver, A)
            % Reciprocal of diagonal matrix
            sz = size(A);
            assert(sz(1) == sz(2), 'Matrix must be square!');
            n = sz(1);
            d = 1./abs(diag(A));
            d(~isfinite(d)) = 1;
            I = (1:n)';
            M = sparse(I, I, d, n, n);
        end
        
        function initialGuess = getInitialGuess(solver, problem)
            if isempty(solver.initialGuessFun)
                initialGuess = {[]};
                return
            end
            initialGuess = solver.initialGuessFun(problem);
        end
        
        function [d, shortname] = getDescription(solver)
            % Get the description and a short name used for display
            % purposes.
            shortname = 'Virtual base';
            shortname = [shortname, solver.id];

            d = 'Virtual base class linear solver. Not for direct use.';
        end
        
        function disp(solver)
            [d, sn] = solver.getDescription();
            s1 = sprintf('  %s linear solver of class %s', sn, class(solver));
            fprintf('%s\n  %s\n  %s\n  ->', s1, repmat('-', 1, numel(s1)-2), d);
            builtin('disp', solver);
        end
    end
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
