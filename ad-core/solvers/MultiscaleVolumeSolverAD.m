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
       
       function [result, report] = solveLinearSystem(solver, A, b)
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
    restrict = sparse(j, i, v, n, m);
end
