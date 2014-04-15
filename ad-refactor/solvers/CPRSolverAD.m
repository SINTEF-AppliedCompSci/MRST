classdef CPRSolverAD < linearSolverAD
    properties
        relativeTolerance
        pressureScaling
        ellipticSolver
    end
    methods
        function solver = CPRSolverAD(varargin)
            solver.ellipticSolver = mldivideSolverAD();
            solver.relativeTolerance = 1e-3;
            solver.pressureScaling = 1/(200*barsa);
%             solver.pressureScaling = 1;

            
            disp('OK!')
        end
        
        function result = solveLinearSystem(solver, A, b) %#ok
            error('Not supported')
        end
        
        function [dx, result] = solveLinearProblem(solver, problem)

            
            % Solve a linearized problem
            isPressure = problem.indexOfPrimaryVariable('pressure');
            pressureIndex = find(isPressure);
            
            
            edd = 1e-2;
            
            % Eliminate the non-cell variables first
            isCell = problem.indexOfType('cell');
            cellIndex = find(isCell);
            cellEqNo = numel(cellIndex);
            
            % Find number of "cell" like variables
            n = problem.getEquationVarNum(cellIndex(1));
            nP = numel(problem);
            
            % Eliminate non-cell variables (well equations etc)
            notCellIndex = find(~isCell);
            
            eliminated = cell(numel(notCellIndex), 1);
            elimNames = problem.equationNames(notCellIndex);
            
            for i = 1:numel(notCellIndex)
                [problem, eliminated{i}] = problem.eliminateVariable(elimNames{i});
            end
            
            
            
            eqs = problem.equations;
            diagMod = zeros(n, cellEqNo);
            for i = 1:cellEqNo
                % Find the derivative of current cell block w.r.t the
                % pressure
                pressureJacobi = eqs{cellIndex(i)}.jac{isPressure};
                pressureDiag  = diag(pressureJacobi);
                
                sod = sum(abs(pressureJacobi), 2) - abs(pressureDiag);
                % Find "bad" equations, i.e. equations where the current
                % pressure index does not give a good elliptic jacobian
                diagMod(pressureDiag./sod > edd, i) = pressureIndex;
            end
            
            diagMod(all(diagMod == 0, 2), 1) = 1;
            bad = diagMod(:, 1) == 0;
            
            if any(bad)
                % Switch equations for non-elliptic jacobian components with
                % some other equation that has an elliptic pressure jacobian
                [r, c] = find(diagMod(bad, :));
                firstIndex = accumarray(r,c,[sum(bad), 1], @min);
                diagMod(bad, 1) = firstIndex;
                diagMod = diagMod(:, 1);
                
                
                for i = unique(diagMod)'
                    if i == pressureIndex
                        continue
                    end
                    % Standard switching with some overloaded indexing in
                    % the ad objects
                    isCurrent = diagMod == i;
                    tmp = eqs{pressureIndex}(isCurrent);
                    eqs{pressureIndex}(isCurrent) = eqs{i}(isCurrent);
                    eqs{i}(isCurrent) = tmp;
                end
            end
            
            ok = diagMod(:, 1) == pressureIndex;
            ok = ok | true;
            for i = 1:cellEqNo
                % Add together all equations to create the pressure-like
                % equations that captures all jacobians
                if i == pressureIndex
                    continue
                end
                eqs{pressureIndex}(ok) = eqs{pressureIndex}(ok) + eqs{i}(ok);
            end
            problem.equations = eqs;
            
            % Set up storage for all variables, including those we
            % eliminated previously
            dx = cell(nP, 1);
            
            % Recover non-cell variables
            recovered = false(nP, 1);
            recovered(cellIndex) = true;
            
            % Solve cell equations
                        
%             refSolver = mldivideSolverAD();
%             dxCell = refSolver.solveLinearProblem(problem);
            
            if solver.pressureScaling ~= 1
                for i = 1:numel(problem)
                    problem.equations{i}.jac{pressureIndex} =...
                    problem.equations{i}.jac{pressureIndex}./solver.pressureScaling;
                end
            end

            % ILU0 preconditioner
            [A, b] = problem.getLinearSystem();
            
            [L, U] = ilu(A, struct('type', 'nofill'));

            Ap = problem.equations{pressureIndex}.jac{pressureIndex};
            % We have only cell variables present, and these will have
            % offsets of cellnum long each
            pInx = (1:n).' + n*(pressureIndex-1);
            % solver
            ellipSolve = @(A, b) solver.ellipticSolver.solveLinearSystem(A, b);

            prec = @(r) applyTwoStagePreconditioner(r, A, L, U, Ap, pInx, ellipSolve);
            [cprSol, fl, relres, its] = gmres(A, b, [], solver.relativeTolerance, 40, prec);
            if solver.pressureScaling ~= 1;
                cprSol(pInx) = cprSol(pInx) ./ solver.pressureScaling;
            end
            dispif(mrstVerbose(), 'GMRES converged in %d iterations\n', its(2));
            dxCell = solver.storeIncrements(problem, cprSol);
            
            dx(recovered) = dxCell;
            
            %
            for i = numel(eliminated):-1:1
                pos = notCellIndex(i);
                % OBS: cellEqNo + 1 is a hack WHICH ASSUMES THAT ALL
                % ELIMINATED VARIABLES COMES AFTER CELL VARIABLES
                dVal = recoverVars(eliminated{i}, cellEqNo + 1, dx(recovered));
                dx{pos} = dVal;
                recovered(pos) = true;
            end
            
            if nargout > 1
                result = vertcat(dx{:});
            end
        end
        
    end
end

function x = applyTwoStagePreconditioner(r, A, L, U, Ap, pInx, ellipticSolver)
   x = zeros(size(r));
   x(pInx) = ellipticSolver(Ap, r(pInx));

   r = r - A*x;
   x = x + U\(L\r);
end
