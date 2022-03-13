classdef MultiscaleVolumeSolverAD < LinearSolverAD
    % Multiscale linear solver
   properties
       prolongationType
       controlVolumeRestriction
       updateBasis
       resetBasis
       basis
       localSolver
       coarsegrid
       setupTime
       useMEX
       mexGrid
       basisIterations
       basisTolerance
       updateCounter
       updateInterval
       
       getSmoother
       useGMRES
   end
   
   properties (Access = private)
       smoother_fn = [];
   end
   methods
       function solver = MultiscaleVolumeSolverAD(coarsegrid, varargin)
           solver = solver@LinearSolverAD();
           
           Nf = coarsegrid.parent.cells.num;
           Nc = coarsegrid.cells.num;
           dim = coarsegrid.parent.griddim;
           
           % Default options
           solver.prolongationType = 'msrsb';
           solver.controlVolumeRestriction = true;
           solver.updateBasis = false;
           solver.maxIterations = 0;
           solver.getSmoother = [];
           solver.useGMRES = false;
           solver.basisTolerance = 5e-3;
           solver.useMEX = true;
           solver.mexGrid = [];
           solver.resetBasis = false;
           solver.updateInterval = 1;
           solver.basisIterations = ceil(50*(Nf/Nc).^(1/dim));
           solver.keepNumber = coarsegrid.parent.cells.num;
           
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
           if isempty(solver.smoother_fn)
               fn = solver.getSmoother;
           else
               fn = solver.smoother_fn;
           end
           [x, report] = solveMultiscaleIteratively(A, b, [], ...
                                                          solver.basis, ...
                                                          fn, ...
                                                          solver.tolerance,...
                                                          solver.maxIterations, ...
                                                          @(A, b) mldivide(A, b), ...
                                                          solver.useGMRES, ...
                                                          solver.verbose);
       end
              
       function solver = setupSolver(solver, A, b, varargin)
           % Run setup on a solver for a given system
           solver = solver.createBasis(A);
           solver.smoother_fn = solver.getSmoother(A, b);
       end
       
       function  solver = cleanupSolver(solver, A, b, varargin) %#ok 
           % Clean up solver after use (if needed)
           solver.smoother_fn = [];
       end
       
       function [dx, result, report] = solveLinearProblem(solver, problem0, model)
           % Solve a linearized problem
           timer = tic();
           problem = problem0.assembleSystem();
           [A, b] = problem.getLinearSystem();
           [A, b, lsys] = solver.reduceLinearSystem(A, b);
           t_prepare = toc(timer);
           if isempty(solver.basis)
               solver = solver.createBasis(A);
           else
               if problem.iterationNo == 1
                   if isempty(solver.updateCounter)
                       solver.updateCounter = 0;
                   end
                   if mod(solver.updateCounter, solver.updateInterval) == 0
                       if solver.resetBasis
                          solver.basis = [];
                          solver = solver.createBasis(A);
                       elseif solver.updateBasis
                          solver = solver.createBasis(A);
                       end
                       solver.updateCounter = 1;
                   else
                       solver.updateCounter = solver.updateCounter + 1;
                   end
               end
           end
           t_basis = toc(timer) - t_prepare;
           [result, report] = solver.solveLinearSystem(A, b);
           t_solve = toc(timer) - t_prepare - t_basis;
           result = solver.recoverLinearSystem(result, lsys);
           
           [result, report] = problem.processResultAfterSolve(result, report);
           report.SolverTime = toc(timer);
           report.LinearSolutionTime = t_solve;
           report.BasisTime = t_basis;
           report.PreparationTime = t_prepare;
           report.PostProcessTime = report.SolverTime - t_solve - t_prepare;
           dx = solver.storeIncrements(problem, result);
       end

       function solver = createBasis(solver, A)
           if isempty(solver.basis) || solver.updateBasis
               if solver.verbose
                   fprintf('Constructing multiscale basis of type %s', solver.prolongationType)
                   if solver.controlVolumeRestriction
                       fprintf(' with control volume (CV) restriction.\n');
                   else
                       fprintf(' with Galerkin (FE) restriction.\n');
                   end
               end
               timer = tic();
               [ii, jj, vv] = find(A);
               % Pick out the diagonal
               d = zeros(size(A, 1), 1);
               isDiag = ii == jj;
               d(ii(isDiag)) = vv(isDiag);
               % Take the positive off-diagonal entries, following the sign
               % convention of the diagonal
               keep = ii ~= jj & sign(d(ii)) == -sign(vv);
               n = size(A, 1);
               % Strip down the sparse structure to the kept elements
               ii = ii(keep);
               jj = jj(keep);
               vv = vv(keep);
               % Take the row sum to be the diagonal (with opposite sign)
               dd = accumarray(ii, vv, [size(A, 1), 1]);
               dd(dd == 0) = mean(dd);
               % Build sparse matrix format
               ii = [ii; (1:n)'];
               jj = [jj; (1:n)'];
               vv = [vv; -dd];
               % This is the regularized matrix used to get
               % partition-of-unity basis functions. Note that there is an
               % alternative reguarlization in getMultiscaleBasis, which we
               % explicitly disable here since it is already fixed.
               A = sparse(ii, jj, vv, n, n);               
               [solver.basis, solver.coarsegrid] =...
                   getMultiscaleBasis(solver.coarsegrid, A, 'useMEX', solver.useMEX, ...
                                                            'mexGrid', solver.mexGrid, ...
                                                            'basis',  solver.basis, ...
                                                            'iterations', solver.basisIterations, ...
                                                            'tolerance', solver.basisTolerance, ...
                                                            'type',   solver.prolongationType, ...
                                                            'useControlVolume', solver.controlVolumeRestriction, ...
                                                            'regularizeSys', false);
               solver.setupTime = toc(timer);
               dispif(solver.verbose, 'Basis constructed in %s.\n', formatTimeRange(solver.setupTime));
           end
       end
       
        function [d, sn] = getDescription(solver)
            sn = solver.prolongationType;
            sn = [sn, solver.id];
            if solver.controlVolumeRestriction
                v = 'volume';
            else
                v = 'element';
            end
            d = sprintf('Multiscale finite-%s solver with %s basis functions', v, solver.prolongationType);
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
