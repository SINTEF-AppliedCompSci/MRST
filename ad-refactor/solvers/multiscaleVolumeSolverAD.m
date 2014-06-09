classdef multiscaleVolumeSolverAD < linearSolverAD
    % Base class for a nonlinear solver
   properties
       prolongationType
       prolongationOperator
       restrictionType
       restrictionOperator
       useGalerkinRestriction
       basis
       localSolver
       coarsegrid
   end
   methods
       function solver = multiscaleVolumeSolverAD(coarsegrid, varargin)
           opt = struct('prolongationType', 'smoothed', ...
                        'restrictionType',  'controlVolume', ...
                        'galerkinRestriction', true);
           opt = merge_options(opt, varargin{:});
           
           solver.prolongationType = opt.prolongationType;
           solver.restrictionType  = opt.restrictionType;
           solver.useGalerkinRestriction = opt.galerkinRestriction;
           solver.coarsegrid = coarsegrid;
       end
       
       function result = solveLinearSystem(solver, A, b)
           I = solver.prolongationOperator;
           
           if solver.useGalerkinRestriction
               R = I';
           else
               R = solver.restrictionOperator;
           end
           coarseValues = (R*A*I)\(R*b);
           result = I*coarseValues;
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