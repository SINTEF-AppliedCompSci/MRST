classdef CPRSolverAD < linearSolverAD
    properties
        relativeTolerance
        pressureScaling
        ellipticSolver
    end
    methods
        function solver = CPRSolverAD(varargin)
            solver = solver@linearSolverAD();
            
            % Default options
            solver.ellipticSolver = [];
            solver.relativeTolerance = 1e-3;
            solver.pressureScaling = 1/(200*barsa);
            
            solver = merge_options(solver, varargin{:});
            
            if isempty(solver.ellipticSolver)
                solver.ellipticSolver = mldivideSolverAD();
            else
                assert(isa(solver.ellipticSolver, 'linearSolverAD'));
            end
        end
        
        function [result, report] = solveLinearSystem(solver, A, b) %#ok
            error('Not supported directly - this is a preconditioner')
        end
        
        function [dx, result, report] = solveLinearProblem(solver, problem, model)
            % Solve a linearized problem using a constrained pressure
            % residual preconditioner
            
            isPressure = problem.indexOfPrimaryVariable('pressure');
            pressureIndex = find(isPressure);
            
            
            edd = 1e-2;
            
            % Eliminate the non-cell variables first
            isCell = problem.indexOfType('cell');
            cellIndex = find(isCell);
            cellEqNo = numel(cellIndex);
            
            % Find number of "cell" like variables
            nCell = problem.getEquationVarNum(cellIndex(1));
            nP = numel(problem);
            
            % Eliminate non-cell variables (well equations etc)
            notCellIndex = find(~isCell);
            
            eliminated = cell(numel(notCellIndex), 1);
            elimNames = problem.equationNames(notCellIndex);
            
            for i = 1:numel(notCellIndex)
                [problem, eliminated{i}] = problem.eliminateVariable(elimNames{i});
            end

            
            eqs = problem.equations;
                        
            % Get and apply scaling
            scale = getScaling(problem, model);
            for eqn = 1:cellEqNo
                eqs{eqn} = eqs{eqn}.*scale(eqn);
            end
            
            
            isElliptic = false(nCell, cellEqNo);
            for i = 1:cellEqNo
                % Find the derivative of current cell block w.r.t the
                % pressure
                pressureJacobi = eqs{i}.jac{isPressure};
                pressureDiag  = diag(pressureJacobi);
                
                sod = sum(abs(pressureJacobi), 2) - abs(pressureDiag);
                % Find "bad" equations, i.e. equations where the current
                % pressure index does not give a good elliptic jacobian
                isElliptic(pressureDiag./sod > edd, i) = true;
            end
            
            isElliptic(all(isElliptic == 0, 2), pressureIndex) = true;
            
            bad = ~isElliptic(:, pressureIndex);

            if any(bad)
                % Switch equations for non-elliptic jacobian components with
                % some other equation that has an elliptic pressure jacobian
                [r, c] = find(isElliptic(bad, :));
                replacementInd = accumarray(r,c,[sum(bad), 1], @min);
                
                
                for i = unique(replacementInd)'
                    if i == pressureIndex
                        continue
                    end
                    % Standard switching with some overloaded indexing in
                    % the ad objects
                    isCurrent = replacementInd == i;
                    tmp = eqs{pressureIndex}(isCurrent);
                    eqs{pressureIndex}(isCurrent) = eqs{i}(isCurrent);
                    eqs{i}(isCurrent) = tmp;
                    
                    isElliptic(isCurrent, pressureIndex) = true;
                    isElliptic(isCurrent, i) = false;
                end
            end
            
            for i = 1:cellEqNo
                % Add together all the equations to get a "pressure
                % equation" where we know that all submatrices should be as
                % close as possible to M-matrices (because of switching,
                % which does not actually alter the solution)
                if i == pressureIndex
                    continue
                end
                ok = isElliptic(:, i);
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
            if solver.pressureScaling ~= 1
                for i = 1:numel(problem)
                    problem.equations{i}.jac{pressureIndex} =...
                    problem.equations{i}.jac{pressureIndex}./solver.pressureScaling;
                end
            end

            % Set up linear system
            [A, b] = problem.getLinearSystem();
            
            % ILU0 preconditioner for the non-elliptic part
            [L, U] = ilu(A, struct('type', 'nofill'));

            Ap = -problem.equations{pressureIndex}.jac{pressureIndex};
            % We have only cell variables present, and these will have
            % offsets of cellnum long each
            pInx = (1:nCell).' + nCell*(pressureIndex-1);
            
            
            % Set up elliptic solver
            solver.ellipticSolver = solver.ellipticSolver.setupSolver(Ap, b(pInx));
            ellipSolve = @(b) solver.ellipticSolver.solveLinearSystem(Ap, b);

            prec = @(r) applyTwoStagePreconditioner(r, A, L, U, pInx, ellipSolve);
            [cprSol, fl, relres, its, resvec] = gmres(A, b, [], solver.relativeTolerance, 40, prec);
            
            % Undo pressure scaling
            if solver.pressureScaling ~= 1;
                cprSol(pInx) = cprSol(pInx) ./ solver.pressureScaling;
            end
            % Clean up elliptic solver
            solver.ellipticSolver = solver.ellipticSolver.cleanupSolver(Ap, b(pInx));
            
            dispif(mrstVerbose(), 'GMRES converged in %d iterations\n', its(2));
            dxCell = solver.storeIncrements(problem, cprSol);
            
            % Put the recovered variables into place
            dx(recovered) = dxCell;
            
            assert(all(diff(cellIndex) == 1), 'This solver currently assumes that the cell variables comes first!')
            for i = numel(eliminated):-1:1
                pos = notCellIndex(i);
                dVal = recoverVars(eliminated{i}, cellEqNo + 1, dx(recovered));
                dx{pos} = dVal;
                recovered(pos) = true;
            end
            
            if nargout > 1
                result = vertcat(dx{:});
            end
            
            if nargout > 2
                report = struct('IterationsGMRES', its(2), ...
                                'FlagGMRES',       fl, ...
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

function scale = getScaling(problem, model)
    state = problem.state;
    fluid = model.fluid;
    p  = mean(state.pressure);
    
    if model.water
        bW = fluid.bW(p);
    end
    
    scale = ones(numel(problem), 1);
    
    isBO = isa(model, 'threePhaseBlackOilModel');
    
    if isBO && model.disgas
        rs = fluid.rsSat(p);
        bO = fluid.bO(p, rs, true);
    elseif model.oil
        bO = fluid.bO(p);
    end
    if isBO && model.vapoil
        rv = fluid.rvSat(p);
        bG = fluid.bG(p, rv, true);
    elseif model.gas
        bG = fluid.bG(p);
    end

    if model.oil
        scale(problem.indexOfEquationName('oil')) = 1./bO;
    end
    if model.gas
        scale(problem.indexOfEquationName('gas')) = 1./bG;
    end
    if model.water
        scale(problem.indexOfEquationName('water')) = 1./bW;
    end
end