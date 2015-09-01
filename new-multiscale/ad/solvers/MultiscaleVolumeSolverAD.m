classdef MultiscaleVolumeSolverAD < LinearSolverAD
    % Multiscale linear solver
   properties
       prolongationType
       controlVolumeRestriction
       updateBasis
       basis
       localSolver
       coarsegrid
       setupTime
       useMEX
       mexGrid
       basisIterations
       
       getSmoother
       useGMRES
   end
   methods
       function solver = MultiscaleVolumeSolverAD(coarsegrid, varargin)
           solver = solver@LinearSolverAD();
           
           Nf = coarsegrid.parent.cells.num;
           Nc = coarsegrid.cells.num;
           dim = coarsegrid.parent.griddim;
           
           % Default options
           solver.prolongationType = 'smoothed';
           solver.controlVolumeRestriction = true;
           solver.updateBasis = false;
           solver.maxIterations = 0;
           solver.getSmoother = [];
           solver.useGMRES = false;
           solver.useMEX = true;
           solver.mexGrid = [];
           solver.basisIterations = ceil(50*(Nf/Nc).^(1/dim));
           
           solver = merge_options(solver, varargin{:});
           
           solver.setupTime = 0;
           solver.coarsegrid = coarsegrid;
       end
       
       function [x, report] = solveLinearSystem(solver, A, b)
           CG = solver.coarsegrid;
           nc = CG.parent.cells.num;
           if isempty(solver.basis)
               solver.setupSolver(A(1:nc, 1:nc), b(1:nc));
           end
           
           [x, report] = solveMultiscaleIteratively(A, b, solver.basis, ...
                                                          solver.getSmoother, ...
                                                          solver.tolerance,...
                                                          solver.maxIterations, ...
                                                          @(A, b) mldivide(A, b), ...
                                                          solver.useGMRES, ...
                                                          solver.verbose);
       end
              
       function solver = setupSolver(solver, A, b, varargin) %#ok 
           % Run setup on a solver for a given system
           solver = solver.createBasis(A);
       end
       
       function  solver = cleanupSolver(solver, A, b, varargin) %#ok 
           % Clean up solver after use (if needed)
       end
       
       function [dx, result, report] = solveLinearProblem(solver, problem0, model)
           % Solve a linearized problem
           [problem, eliminated] = problem0.reduceToSingleVariableType('cell');
           
           problem = problem.assembleSystem();
           
           timer = tic();
           [result, report] = solver.solveLinearSystem(problem.A, problem.b); 
           report.SolverTime = toc(timer);
           
           dxCell = solver.storeIncrements(problem, result);
           dx = problem.recoverFromSingleVariableType(problem0, dxCell, eliminated);
       end

       function solver = createBasis(solver, A)
           if isempty(solver.basis)
               if solver.verbose
                   fprintf('Constructing multiscale basis of type %s', solver.prolongationType)
                   if solver.controlVolumeRestriction
                       fprintf(' with control volume (CV) restriction.\n');
                   else
                       fprintf(' with Galerkin (FE) restriction.\n');
                   end
               end
               timer = tic();
               [solver.basis, solver.coarsegrid] =...
                   getMultiscaleBasis(solver.coarsegrid, A, 'useMEX', solver.useMEX, ...
                                                            'mexGrid', solver.mexGrid, ...
                                                            'iterations', solver.basisIterations, ...
                                                            'type',   solver.prolongationType, ...
                                                            'useControlVolume', solver.controlVolumeRestriction, ...
                                                            'regularizeSys', true);
               solver.setupTime = toc(timer);
               dispif(solver.verbose, 'Basis constructed in %s.\n', formatTimeRange(solver.setupTime));
           elseif solver.updateBasis
               % THIS PART SHOULD BE IMPROVED TO SUPPORT OTHER SOLVERS!
               A_basis = A - diag(sum(A, 2));
               timer = tic();
               switch solver.prolongationType
                   case 'smoothed'
                       % Just redo iterations until we are below increment
                       % tolerance again
                       I = iteratedJacobiBasis(A_basis, solver.coarsegrid,...
                                            'useConstant', true, ...
                                            'interpolator', solver.prolongationOperator);
                       solver.basis.I = I;
                   otherwise
                       error('Error!')
               end
               solver.setupTime = solver.setupTime + toc(timer);
               solver.prolongationOperator = I;
           end
       end
   end
end
