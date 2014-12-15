classdef CPRSolverAD < LinearSolverAD
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
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   See class properties.
%
%
% SEE ALSO:
%   BackslashSolverAD, LinearSolverAD, LinearizedProblem

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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
        % Relative tolerance for elliptic solver
        relativeTolerance
        % Scaling factor applied to pressure equations
        pressureScaling
        % LinearSolverAD subclass suitable for the elliptic submatrix.
        ellipticSolver
        % Diagonal tolerance in [0,1].
        diagonalTol
    end
    methods
        function solver = CPRSolverAD(varargin)
            solver = solver@LinearSolverAD();
            
            % Default options
            solver.ellipticSolver = [];
            solver.relativeTolerance = 1e-2;
            solver.pressureScaling = 1/(200*barsa);
            solver.diagonalTol = 1e-2;
            
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
            
            isPressure = problem0.indexOfPrimaryVariable('pressure');
            pressureIndex = find(isPressure);
            
            [problem, eliminated] = problem0.reduceToSingleVariableType('cell');
            isCurrent = problem.indexOfType('cell');
            cellIndex = find(isCurrent);
            
            cellEqNo = numel(cellIndex);
            nCell = problem.getEquationVarNum(1);
            
            % Get and apply scaling
            eqs = problem.equations;
            scale = getScaling(problem, model);
            for eqn = 1:cellEqNo
                eqs{eqn} = eqs{eqn}.*scale(eqn);
            end
            problem.equations = eqs;
            
            isElliptic = false(nCell, cellEqNo);
            for i = 1:cellEqNo
                % Find the derivative of current cell block w.r.t the
                % pressure
                pressureJacobi = eqs{i}.jac{isPressure};
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
                %[r, c] = find(isElliptic(bad, :));
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
            for i = 1:cellEqNo
                % Add together all the equations to get a "pressure
                % equation" where we know that all submatrices should be as
                % close as possible to M-matrices (because of switching,
                % which does not actually alter the solution)
                ok = isElliptic(:, i);
                if i == pressureIndex
                    continue
                end
                problem.equations{pressureIndex} = ...
                    problem.equations{pressureIndex} + ok.*eqs{i};
            end
            

            
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
            try
                [cprSol, fl, relres, its, resvec] = gmres(A, b, [], solver.relativeTolerance,...
                                                    min(solver.maxIterations, size(A, 1)), prec);
            catch exception
                % Ensure external memory etc is deallocated properly if the
                % solve failed in some spectacular manner
                solver.ellipticSolver = solver.ellipticSolver.cleanupSolver(Ap, b(pInx));
                rethrow(exception)
            end
            % Clean up elliptic solver
            solver.ellipticSolver = solver.ellipticSolver.cleanupSolver(Ap, b(pInx));

            % Undo pressure scaling
            if solver.pressureScaling ~= 1;
                cprSol(pInx) = cprSol(pInx) ./ solver.pressureScaling;
            end
            
            
            dxCell = solver.storeIncrements(problem, cprSol);
            dx = problem.recoverFromSingleVariableType(problem0, dxCell, eliminated);
            
            %dispif(solver.verbose, 'GMRES converged in %d iterations\n', its(2));
            % Recover stuff
            solvetime = toc(timer);
            
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

function scale = getScaling(problem, model)
    state = problem.state;
    fluid = model.fluid;
    p  = mean(state.pressure);
    
    if model.water
        bW = fluid.bW(p);
    end
    
    scale = ones(numel(problem), 1);
    
    isBO = isa(model, 'ThreePhaseBlackOilModel');
    
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
