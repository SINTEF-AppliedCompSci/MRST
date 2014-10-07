classdef MultiscaleVolumeSolverAD < LinearSolverAD
    % Multiscale linear solver
   properties
       prolongationType
       prolongationOperator
       restrictionType
       restrictionOperator
       useGalerkinRestriction
       updateBasis
       basis
       localSolver
       coarsegrid
       setupTime
   end
   methods
       function solver = MultiscaleVolumeSolverAD(coarsegrid, varargin)
           solver = solver@LinearSolverAD();
           
           % Default options
           solver.prolongationType = 'smoothed';
           solver.restrictionType  = 'controlVolume';
           solver.useGalerkinRestriction = true;
           solver.updateBasis = false;
           
           solver = merge_options(solver, varargin{:});
           
           solver.setupTime = 0;
           solver.coarsegrid = coarsegrid;
       end
       
       function [dx, result, report] = solveLinearProblem(solver, problem, model)
           
            % Eliminate the non-cell variables first
            isCell = problem.indexOfType('cell');
            cellIndex = find(isCell);
            cellEqNo = numel(cellIndex);
            
            % Find number of "cell" like variables
            nCell = problem.getEquationVarNum(cellIndex(1));
            nP = numel(problem);
                        
            % Eliminate non-cell variables (well equations etc)
            problem = problem.clearSystem();
            
            notCellIndex = find(~isCell);
            
            eliminated = cell(numel(notCellIndex), 1);
            elimNames = problem.equationNames(notCellIndex);
            
            for i = 1:numel(notCellIndex)
                [problem, eliminated{i}] = problem.eliminateVariable(elimNames{i});
            end
            % SOLVE HERE
           % Solve a linearized problem
           problem = problem.assembleSystem();
           
           timer = tic();
           [result, report] = solver.solveLinearSystem(problem.A, problem.b); 
%            tmp = BackslashSolverAD();
%            [result, report] = tmp.solveLinearSystem(problem.A, problem.b); 
           report.SolverTime = toc(timer);

           dxCell = solver.storeIncrements(problem, result);
           
            % Set up storage for all variables, including those we
            % eliminated previously
            dx = cell(nP, 1);
            
            % Recover non-cell variables
            recovered = false(nP, 1);
            recovered(cellIndex) = true;
            
            % Put the recovered variables into place
            dx(recovered) = dxCell;
            
            for i = numel(eliminated):-1:1
                pos = notCellIndex(i);
                dVal = recoverVars(eliminated{i}, cellEqNo + 1, dx(recovered));
                dx{pos} = dVal;
                recovered(pos) = true;
            end
       end
       
       function [result, report] = solveLinearSystem(solver, A, b)
           if isempty(solver.prolongationOperator)
               solver.setupSolver(A, b);
           end
           
           I = solver.prolongationOperator;
           
           if solver.useGalerkinRestriction
               R = I';
           else
               R = solver.restrictionOperator;
           end
           coarseValues = (R*A*I)\(R*b);
           result = I*coarseValues;
           % Nothing to report
           report = struct();
       end
       
       function solver = setupSolver(solver, A, b, varargin) %#ok 
           % Run setup on a solver for a given system
           solver = solver.createBasis(A);
           if ~solver.useGalerkinRestriction
               solver.restrictionOperator =...
                    restrictOperator(solver.coarsegrid.partition);
           end
       end
       
       function  solver = cleanupSolver(solver, A, b, varargin) %#ok 
           % Clean up solver after use (if needed)

       end
       
       function solver = createBasis(solver, A)
           if isempty(solver.basis)
               A_basis = A - diag(sum(A, 2));
               disp('Constructing basis!')
               timer = tic();
               switch solver.prolongationType
                   case 'mstpfa'
                        fb = createFaceBasis(CG, A_basis);
                        I = assembleCoarseOperatorsPartition(CG, fb);
                        solver.basis = fb;
                   case 'smoothed'
                       I = iteratedJacobiBasis(A_basis, solver.coarsegrid,...
                                            'useConstant', true);
                       solver.basis = I;
                   otherwise
                       error('Error!')
               end
               disp('Basis constructed!')
               solver.setupTime = toc(timer);
               solver.prolongationOperator = I;
           elseif solver.updateBasis
               A_basis = A - diag(sum(A, 2));
               timer = tic();
               switch solver.prolongationType
                   case 'smoothed'
                       % Just redo iterations until we are below increment
                       % tolerance again
                       I = iteratedJacobiBasis(A_basis, solver.coarsegrid,...
                                            'useConstant', true, ...
                                            'interpolator', solver.prolongationOperator);
                       solver.basis = I;
                   otherwise
                       error('Error!')
               end
               solver.setupTime = solver.setupTime + toc(timer);
               solver.prolongationOperator = I;
           end
       end
   end
end


function restrict = restrictOperator(partition)
    n = numel(partition);
    m = max(partition);
    i = (1:n) .';
    j = partition(i);
    v = ones(n, 1);
    restrict = sparse(i, j, v, n, m)';
end
