function state = solveMSFV_TPFA_Incomp_DFM(state, G, CG, T, fluid, varargin)
%SOLVEMSFV_TPFA_INCOMP_DFM Modified to allow for hybrid cells
% finite volume method.
%
% THIS FILE IS MODIFIED FROM THE ORIGINAL solveMSFV_TPFA_Incomp.m
% to allow for hybrid cells and cell to cell connections
%
% The change are:
%   - optional input c2cTrans is added to allow for c2c connections
%     This is used as an input in incompTPFA_DFM
%   - incompTPFA_DFM is called instead of incompTPFA to allow for
%     hybrid cells
%   - a bug in the flux reconstruction is fixed
%
% Portions Copyright 2013 IRIS AS
%
%
% SYNOPSIS:
% state = solveMSFV_TPFA_Incomp(state, G, CG, T, fluid)
% state = solveMSFV_TPFA_Incomp(state, G, CG, T, fluid, 'pn1', pv1, ...)
%
% DESCRIPTION:
% This function uses a operator formulation of the multiscale finite volume
% method to solve a single phase flow problem. The method works on fully
% unstructured grids, as long as the coarse grids are successfully created.
%
% REQUIRED PARAMTERS:
%   state   - reservoir and well solution structure as defined by
%             'initResSol', 'initWellSol' or a previous call to 'incompTPFA'.
%
%   G, T   - Grid and half-transmissibilities as computed by the function
%            'computeTrans'.
%
%   fluid  - Fluid object as defined by function 'initSimpleFluid'.
%
% OPTIONAL PARAMETERS:
%   wells  - Well structure as defined by function 'addWell'.  May be empty
%            (i.e., W = struct([])) which is interpreted as a model without
%            any wells.
%
%   bc     - Boundary condition structure as defined by function 'addBC'.
%            This structure accounts for all external boundary conditions to
%            the reservoir flow.  May be empty (i.e., bc = struct([])) which
%            is interpreted as all external no-flow (homogeneous Neumann)
%            conditions.
%
%   src    - Explicit source contributions as defined by function
%            'addSource'.  May be empty (i.e., src = struct([])) which is
%            interpreted as a reservoir model without explicit sources.
%
%   LinSolve     - Handle to linear system solver software which will be
%                  used to solve the linear equation sets.
%                  Assumed to support the syntax
%
%                        x = LinSolve(A, b)
%
%                  in order to solve a system Ax=b of linear equations.
%                  Should be capable of solving for a matrix rhs.
%                  Default value: LinSolve = @mldivide (backslash).
%
%   Dual        - A dual grid. Will be created if not already defined.
%
%   Verbose     - Controls the amount of output. Default value from
%                 mrstVerbose global
%
%   Reconstruct - Determines if a conservative flux field will be
%                 constructed after initial pressure solution.
%
%   Discretization - Will select the type of stencil used for approximating
%                    the pressure equations. Currently only TPFA.
%
%   Iterations     - The number of iterations performed. With subiterations
%                    and a smoother, the total will be subiter*iter.
%
%   Subiterations  - When using smoothers, this gives the number of
%                    smoothing steps.
%
%   Tolerance      - Termination criterion for the iterative variants.
%
%   CoarseDims     - For logical partitioning scheme, the coarse dimensions
%                    must be specified.
%
%   Restart        - The number of steps before GMRES restarts.
%
%   Iterator       - GMRES, DAS or DMS. Selects the iteration type.
%
%   Omega          - Relaxation paramter for the smoothers
%
%   Scheme         - The scheme used for partitioning the dual grid. Only
%                    used when no dual grid is supplied.
%
%   DoSolve        - Will the linear systems be solved or just outputted?
%
%   Speedup        - Alternative formulation of the method leads to great
%                    speed improvements in 3D, but may have different error
%                    values and is in general more sensitive to grid errors
%
%   Update         - If a state from a previous multiscale iteration is
%                    provided, should operators be regenerated?
%
%   DynamicUpdate  - Update pressure basis functions adaptively. Requires
%                    dfs search implemented as components
% RETURNS:
%   state - Update reservoir and well solution structure with new values
%           for the fields:
%              - pressure -- Pressure values for all cells in the
%                            discretised reservoir model, 'G'.
%              - pressure_reconstructed -- Pressure values which gives
%                                          conservative flux
%              - pressurecoarse -- Coarse pressure used for solution
%              - facePressure --
%                            Pressure values for all interfaces in the
%                            discretised reservoir model, 'G'.
%              - time     -- Timing values for benchmarking
%              - flux     -- Flux across global interfaces corresponding to
%                            the rows of 'G.faces.neighbors'.
%              - M_ee,A_ii-- Relevant matrices for debugging
%              - DG       -- Dual grid used
%              - A        -- System matrix.  Only returned if specifically
%                            requested by setting option 'MatrixOutput'.
%
%              - wellSol  -- Well solution structure array, one element for
%                            each well in the model, with new values for
%                            the fields:
%                              - flux     -- Perforation fluxes through all
%                                            perforations for corresponding
%                                            well.  The fluxes are
%                                            interpreted as injection
%                                            fluxes, meaning positive
%                                            values correspond to injection
%                                            into reservoir while negative
%                                            values mean
%                                            production/extraction out of
%                                            reservoir.
%                              - pressure -- Well bottom-hole pressure.
%               - c2cTrans  -- cell to cell transmissibilities see
%                              computeHybridTrans.m, used as input in
%                              incompTPFA_DFM
%
% NOTE:
%   This solver is based on the core MRST function incompTPFA.
%
% SEE ALSO:
%   `incompTPFA`.

%{
Copyright 2009, 2010, 2011, 2012 SINTEF ICT, Applied Mathematics.

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


    %% Merge defaults with specified arguments
    opt = struct('bc', [], 'src', [], 'wells', [], ...
                'LinSolve',         @mldivide,        ...
                'Dual',             [],...
                'Verbose',          mrstVerbose, ...
                'Reconstruct',      false,...
                'Discretization',   'tpfa', ...
                'Iterations',       0, ...
                'Subiterations',       0, ...
                'Tolerance',        1e-9,...
                'CoarseDims',       [], ...
                'Restart',          [], ...
                'Iterator',         'gmres',...
                'Smoother',         'DMS',...
                'Omega',            5, ...
                'Scheme',           'plane',      ...
                'DoSolve',          true, ...
                'Update',           false, ...
                'DynamicUpdate',    false,...
                'UpdateThreshold',  0.1,...
                'A',                [],...
                'A_res',             [],...
                'rhs',              [],...
                'rhs_res',           [],...
                'Speedup',          true, ...
                'c2cTrans',         [] ...
                );

    opt = merge_options(opt, varargin{:});
    if CG.griddim == 2
        % Not valid in 2D
        opt.Speedup = false;
    end
    if opt.Speedup
        dispif(opt.Verbose, 'Decoupling edge systems...')
    else
        dispif(opt.Verbose, 'Doing original formulation, will not decouple edge systems...');
    end
    %% Assemble a TFPA matrix system stored in state
    if strcmp(opt.Discretization, 'tpfa')
        s = incompTPFA_DFM(state, G, T, fluid, 'c2cTrans',opt.c2cTrans, 'wells', opt.wells, 'src',opt.src, 'bc', opt.bc, 'Verbose', false, 'MatrixOutput', true, 'onlyMatrixOutput' ,false, 'LinSolve', @(A,b) zeros(size(b)) );
    else
        %Doesn't work yet - need to handle boundary conditions by adding
        %them to the permutation matrix
        s = incompMPFA(state, G, T, fluid, 'wells', opt.wells, 'Verbose', false, 'MatrixOutput', true, 'bc', opt.bc);
    end
    if ~isempty(opt.A)
        state.A = opt.A;
    else
        state.A = s.A;
    end
    if ~isempty(opt.rhs)
        state.rhs = opt.rhs;
    else
        state.rhs = s.rhs;
    end

    clear s;


    if ~isfield(CG.faces, 'centroids')
        % We need expanded coarse geometry to partition the grid
        CG = coarsenGeometry(CG);
    end

    %% Partition the dual grid
    dispif(opt.Verbose, '\nGenerating dual partition...\n');
    t1 = tic;
    if ~isstruct(opt.Dual)
        % If no dual grid has been created, make one
        if strcmpi(opt.Scheme, 'logical')
            assert(~isempty(opt.CoarseDims),'To use the logical partitioning scheme CoarseDims is mandidatory');
            DG = partitionUIdual(CG, opt.CoarseDims);
        elseif strcmpi(opt.Scheme, 'globalplane')
            DG = dualPartitionGlobal(CG);
        else
            DG = dualPartition(CG, 'Verbose', opt.Verbose, 'Scheme', opt.Scheme, 'StoreEdge', opt.Speedup);
        end
    else
        DG = opt.Dual;
    end
    state.time.dualpartition = toc(t1);
    tocif(opt.Verbose, t1);
    %% Generate permutation matrix
    dispif(opt.Verbose, 'Generating permutation matrix...\n');
    [CG DG] = createPermutationMatrix(state.A, DG, CG, 'wells', opt.wells, 'Speedup', opt.Speedup);

    %% Find the various blocks of the system with notation defined in
    %  "Operator MSFVM formulation with correction function" I. Lunati & S.
    %  Lee 2009
    % nn refers to coarse center nodes, ss the edges (in 3D), ee the
    % edges (2D) and faces (3D) and ii the inner (the rest)

    n_i = length(DG.ii); n_e = n_i + length(DG.ee); n_l = n_e + length(DG.ll); n_n = n_l + length(DG.nn);
    Nf = DG.N;

    % First create the permuted system and rhs
    % Only take equations and unknowns corresponding to actual cells, and
    % not virtual well cells

    A = DG.P*state.A(1:Nf, 1:Nf)*DG.P';
    r = DG.P*state.rhs(1:Nf);

    state.P = DG.P;

    % Internal nodes' influence on internal nodes
    A_ii = A(1:n_i, 1:n_i);
    % Edge nodes' influence on internal nodes
    A_ie = A(1:n_i, (n_i+1):n_e);
    % Internal nodes' influence on internal nodes
    %A_ei = A((n_i+1):n_e, 1:n_i);
    % Edge nodes' influence on edge nodes
    A_ee = A((n_i+1):n_e, (n_i+1):n_e);
    % Edge nodes' influence on center nodes (not in use)
    %A_ne = A((n_e+1):n_n, (n_i+1):n_e);
    % Center nodes' influence on edge nodes
    A_en = A((n_i+1):n_e, (n_l+1):n_n);
    % Internal nodes' influence on internal nodes (not in use and not meaningful for TPFA)
    %A_nn = A((n_s+1):n_n, (n_s+1):n_n);
    if opt.Speedup
        A_ll = A(n_e+1:n_l, n_e+1:n_l);
        A_le = A(n_e+1:n_l, n_i+1:n_e);
        A_li = A(n_e+1:n_l, 1:n_i    );
        A_ln = A(n_e+1:n_l, n_l+1:n_n);
        A_el = A(n_i+1:n_e, n_e+1:n_l);
        M_ll = A_ll + diag(sum(A_le, 2)) + diag(sum(A_li, 2));
        clear A_ll;
        Ns = length(M_ll);
        if isinf(condest(M_ll))
            warning('msfv:dcpl','Bad decoupling of linear edges')
                % Singular values sometimes appear for non-convex cells. This does not
                % affect the rest of the grid and has already been warned about
                warning off
        end
    end

    %% Generate operators and inverse matrices needed
    dispif(opt.Verbose, 'Generating operators...\n');
    t1 = tic;

    X = restrictOperator(CG, DG, Nf);
    M_ee = A_ee + diag(sum(A_ie, 1));

    state.M_ee = M_ee;
    state.A_ii = A_ii;
    if opt.Speedup
       state.M_ll = M_ll;
    end
    Ni = length(A_ii);
    Ne = length(M_ee);
    if opt.Verbose
        fprintf('Condition numbers A_ii: %s, M_ee %s\n', num2str(condest(A_ii)), num2str(condest(M_ee)));
    end

    if ~opt.DoSolve
       return
    end

    %% B represents the pressure basis functions, which can be reused
    if isfield(state, 'msfvm')  && ~(opt.Update || opt.DynamicUpdate)
        dispif(opt.Verbose, 'Reusing existing basis functions...');
        if opt.DynamicUpdate && isfield(state, 'olds')
            if ~isfield(state, 'components')
               % The state contains component listing of all fine cells
               % We need this because the dual grid is implicitly defined
               try
                   state.clusters = find_clusters(A_ii, M_ee);
               catch MatlabException
                   error('Missing matlab_bgl library')
               end
            end
            comp = state.clusters;
            s_ratio = state.s(:,1)./state.olds(:,1);
            e_s = opt.UpdateThreshold;
            targets =  (1/(1+e_s) < s_ratio) & (s_ratio < 1 + e_s);
            if ~any(targets == 1)
                B = state.msfvm.B;
            else
                t_cells = find(targets == 1);
                dispif(opt.Verbose, sprintf('Updating basis functions - %d of %d cells\n', numel(t_cells), G.cells.num));
                t_ii = find(ismember(DG.ii, t_cells));
                inv_A_ii = @(b) mldivide_update(A_ii, b, state.msfvm.B(1:n_i,:), t_ii, comp{1});

                B = updateB(A_ie, A_en, inv_A_ii, @(b) M_ee\b);
                % itererer over t_i, osv
            end
        else
            B = state.msfvm.B;
        end
    else
        if opt.Speedup
            B = formBspeedup(A_ie, A_ii, M_ee, M_ll, A_ln, A_el, opt.LinSolve);
        else
            B = formB(A_ie, A_en, A_ii, M_ee);
        end
        state.msfvm.B = B;
    end
    %% The other operators does not involve caching anything yet
    % Cxr is a ripe target since it does several A_ii inverts, but this is
    % not trivial to implement cleanly


    if opt.Speedup
        Ctimes = @(r) CxrSpeedup(A_ii, M_ee, A_ie, M_ll, A_el, Nf, Ni, Ns, Ne, r, opt.LinSolve);
    else
        Ctimes = @(r) Cxr(A_ii, M_ee, A_ie, Nf, Ni, Ne, r);
    end

    Cr = Ctimes(r);
    dispif(opt.Verbose, 'Ok, operators generated\n');
    state.time.operators = toc(t1);
    tocif(opt.Verbose, t1);

    %% Solving coarse system
    dispif(opt.Verbose, 'Solving coarse system...\n');

    q_n = X*r - X*A*Cr;

    t1 = tic;
    M_nn = X*A*B;
    U_n = mldivide(M_nn, q_n);

    state.M_nn = M_nn;
    state.time.coarsesystem = toc(t1);
    tocif(opt.Verbose, t1);

    %% Interpolate solution and add inn correction functions
    dispif(opt.Verbose, 'Assembling final system...\n');
    t1 = tic;

    U = B*U_n + Cr;

    state.msfvm.Cr = Cr;
    state.time.interpolate = toc(t1);
    tocif(opt.Verbose, t1);

    %% iterations
    if opt.Iterations
        t1 = tic;
        dispif(opt.Verbose, 'Iterations are on\n');
        Nn = CG.cells.num;
        % R: Sparse equivalent of R = [zeros(Nn,Nf-Nn) speye(Nn)];
        R = sparse(1:Nn, Nf - Nn +  (1:Nn),1);
        omega = opt.Omega;
        invG = @(U) Ginv(U, Ctimes, M_nn, Nf, R, X, A, B);

        switch lower(opt.Iterator)
            case 'msfvm'
                % Ordinary MSFVM iterations
                dispif(opt.Verbose, sprintf('Doing MsFV %d iterations with %d sub-smoothing iterations (%s)\n', opt.Iterations, opt.Subiterations, opt.Smoother));
                if opt.Verbose; h = waitbar(0,'Starting iterations...'); end;
                % Create modified systems for the Dirichlet smoothers
                Abar = DG.P_flux*state.A(1:Nf, 1:Nf)*DG.P_flux';
                [D Up] = DU(Abar, CG, Nf);
                rbar = DG.P_flux*(DG.P'*r);
                switch lower(opt.Smoother)
                    case 'dms'
                        smoother = @(res) (D+Up)\res;
                    case 'das'
                        smoother = @(res) D\res;
                    otherwise
                        error('not implemented')
                end
                err = @(Ubar)  norm(invG(r - A*(DG.P*(DG.P_flux'*Ubar))))/norm(invG(rbar));
                last = sum(abs(r - A*U));
                v = 1 ;
                U_new = U;
                while v <= opt.Iterations
                    % Calculate the residual
                    res = r - A*U_new;
                    current = sum(abs(res));
                    if last < current
                        v = v - 1;
                        omega = omega/2;
                        dispif(opt.Verbose, sprintf('Rejecting step, halving omega: %s\n', omega));
                    else
                        U = U_new;
                    end
                    last = current;
                    % Permute to the coarse block ordering via the original
                    % ordering
                    Ubar = DG.P_flux*(DG.P'*U);
                    dispif(opt.Verbose, sprintf('Iteration %d:\n', v));
                    for sub = 1:opt.Subiterations
                        e = err(Ubar);
                        state.residuals((v-1)*opt.Subiterations + sub) = e;
                        dispif(opt.Verbose, sprintf('\tSubstep %d: E = %s\n', sub, e));
                        res = rbar - Abar*Ubar;
                        Ubar = Ubar + (smoother(res));
                    end
                    %Commented because the preconditioned error calculation is really
                    %expensive
%                     if err(Ubar)<opt.Tolerance
%                         break
%                     end
                    % Permute back
                    U = DG.P*(DG.P_flux'*Ubar);
                    % Perform MsFV iteration
                    U_new = U + omega*invG(res);
                    if opt.Verbose;waitbar(v/opt.Iterations,h,sprintf('%d of %d iterations done',(v-1)*opt.Subiterations + sub, opt.Iterations*opt.Subiterations)); end;
                    v = v + 1;
                end
                if opt.Verbose; close(h); end;
            case 'gmres'
                % Generalized Minimized Residual
                dispif(opt.Verbose, sprintf('Doing GMRES with max iterations %d and restart after %d iterations\n', opt.Iterations, opt.Restart));
%                 U = gmres(@(u) A*u,r,opt.Restart,10e-9,opt.Iterations,invG, [], U);
                if isempty(opt.A_res)
                    rr = @(u) A*u;
                else
                    tmpA = DG.P*opt.A_res(1:Nf, 1:Nf)*DG.P';
                    rr = @(u) tmpA*u;
                end

                if ~isempty(opt.rhs_res)
                    r = DG.P*opt.rhs_res(1:Nf);
                end

                t = tic;

                % we use gmres as an inexact solver and do not need it to
                % convergence to the given tolerance.
                warning('off', 'MATLAB:gmres:tooSmallTolerance');
                [U, ~, ~, iter, resvec] = gmres(rr, r ,opt.Restart,opt.Tolerance,opt.Iterations,invG, [], U);
                warning('on', 'MATLAB:gmres:tooSmallTolerance');
                state.t = toc(t);
                state.residuals = resvec./norm(invG(r));
                state.iter = iter(2);
            case 'cg'
                % Generalized Minimzed Residual
                dispif(opt.Verbose, sprintf('Doing CG with max iterations %d\n', opt.Iterations));
                [U, flag, relres, iter, resvec] = pcg(@(u) A*u,r,opt.Tolerance,opt.Iterations,invG, [], U);
                state.residuals = resvec./norm(invG(r));
                state.iter = iter;
            case 'gmres-direct'
                % Generalized Minimzed Residual without preconditioner for
                % comparison purposes
                dispif(opt.Verbose, sprintf('Doing GMRES with max iterations %d and restart after %d iterations (No preconditioner)\n', opt.Iterations, opt.Restart));
%                 U = gmres(@(u) A*u,r,opt.Restart,10e-9,opt.Iterations,invG, [], U);

                [U, flag, relres, iter, resvec] = gmres(@(u) A*u,r,opt.Restart,opt.Tolerance,opt.Iterations,[], []);

                state.residuals = resvec./norm(r);
            otherwise
                error('Unknown solver...')
        end
        state.time.iterations = toc(t1);
        tocif(opt.Verbose, t1);
    end

    p = DG.P'*U;
    state.pressure = p;

    % Dump variables for debugging and plotting
    state.pressurecoarse = U_n;
    state.DG = DG;

    %% Reconstruct a pressure solution with continuous flux from the multiscale solution

    if opt.Reconstruct

        t1 = tic;
        Abar = DG.P_flux*state.A(1:Nf, 1:Nf)*DG.P_flux';
        D = formD(Abar, CG, Nf);

        state.D = D;
        if opt.Verbose
            fprintf('Condition number D: %s\n', num2str(condest(D)));
        end

        % compute the boundary conditions
        q = DG.P_flux*state.rhs(1:Nf)-(Abar - D)*DG.P_flux*state.pressure;

        % solve the localized problems
        % as we have a pure Neumann problem the problem may be singular
        warning('off','MATLAB:nearlySingularMatrix')
        warning('off','MATLAB:singularMatrix')
        sp = mldivide(D,q);
        warning('on','MATLAB:nearlySingularMatrix')
        warning('on','MATLAB:singularMatrix')

        % map back to the orignal indices and update the structure
        p(1:numel(sp)) = DG.P_flux'* sp;
        state.pressure_reconstructed = p(1:G.cells.num);
        state.time.reconstruction = toc(t1);
        tocif(opt.Verbose, t1);
    end

    % Add inn BHP pressures to the solution to reuse existing TPFA code
    tmp = [];
    for i = 1:length(opt.wells)
        w = opt.wells(i);
        if strcmp(w.type,'bhp')
            tmp = [tmp; w.val]; %#ok
        end
    end
    p = vertcat(p, tmp);

    % the original flux is used between the primal cells
    % so we need to compute them as well

    % compute the flux from the reconstructed pressures
    i = all(G.faces.neighbors ~= 0,2);
    state_o = setFlux(G, vertcat(state.pressure,tmp), T, fluid, state, opt);

    % compute the flux from the original pressures
    state = setFlux(G, p, T, fluid, state, opt);

    % change the fluxes
    flux = state.flux;
    flux(CG.faces.fconn) = state_o.flux(CG.faces.fconn);
    state.flux = flux;


    %% check for conservation

    % sum the original fluxes
    sumf_o = accumarray([G.faces.neighbors(i,1);G.faces.neighbors(i,2)],[state_o.flux(i);-state_o.flux(i)]);

    % sum the reconstruced fluxes
    sumf = accumarray([G.faces.neighbors(i,1);G.faces.neighbors(i,2)],[flux(i);-flux(i)]);

    % the fluxes should sum up to the right hand side
    rhs = state.rhs(1:numel(sumf));
    if ~isempty(opt.wells)
        rhs([opt.wells.cells]') = [state.wellSol.flux]';
    end

    % report the improved conservation
    if opt.Verbose
        disp(['Norm of mass conservation error ' num2str(norm(sumf - rhs)) ...
            ' improved from ' num2str(norm(sumf_o - rhs))]);
    end

    % check for coarse conservation
    ic = all(CG.faces.neighbors ~= 0,2);
    cflux = coarsenFlux(CG,flux);
    sumf_c = accumarray([CG.faces.neighbors(ic,1);CG.faces.neighbors(ic,2)],[cflux(ic);-cflux(ic)]);

    % display a warning
    if ( norm(sumf_c - X(:,1:numel(rhs)) * ((DG.P(1:numel(sumf),1:numel(sumf))*rhs))) > 1e-12)
        warning(['Norm of coarse conservation error is larger then ' ...
            num2str(norm(sumf_c - X(:, 1:numel(rhs)) * rhs))]);
    end

    %% Remove virtual entries and save state
    state.pressure = state.pressure(1:G.cells.num);

    %% Save saturation for the next iteration
    state.olds = state.s;
    warning on
    return
end

%% Operator generation
% function M_ee = multiDiagonal(A_ee, A_ie)
%     M_ee = A_ee + diag(sum(A_ie, 1));
% end

function Cr = Cxr(A_ii, M_ee, A_ie, Nf, Ni, Ne, r)
    %multiply C operator with a given vector r optimized
    Cr = zeros(Nf,1);
    Cr(1:Ni) = A_ii\(r(1:Ni) + -A_ie*(M_ee\r((Ni+1):(Ni+Ne))));
    Cr(Ni+1:Ni + Ne) = M_ee\r((Ni+1):(Ni+Ne));
end

function Cr = CxrSpeedup(A_ii, M_ee, A_ie, M_ll, A_el, Nf, Ni, Nl, Ne, r, lsolve)
    % This is the same as Cr, but with an additional interpolation step
    Cr = zeros(Nf,1);
    Cr(1:Ni) = lsolve(A_ii, r(1:Ni) - A_ie*lsolve(M_ee,r((Ni+1):(Ni+Ne)))  + A_ie*lsolve(M_ee,A_el*lsolve(M_ll,r(Ni+Ne+1:Ni+Ne+Nl))));
    Cr(Ni+1:Ni + Ne) = lsolve(M_ee,(r((Ni+1):(Ni+Ne)) -A_el*lsolve(M_ll,r(Ni+Ne+1:Ne+Ni+Nl))));

    Cr(Ni+Ne+1:Ni+Ne+Nl) = lsolve(M_ll,r(Ni+Ne+1:Ni+Ne+Nl));
end

function B = formB(A_ie, A_en, A_ii, M_ee)
    tmp = (M_ee\A_en);
    B = [(A_ii\(A_ie*tmp));...
        -tmp;...
        eye(size(tmp,2))];
    return
end

function B = updateB(A_ie, A_en, inv_A_ii, inv_M_ee)
    tmp = (inv_M_ee(A_en));
    B = [(inv_A_ii(A_ie*tmp));...
        -tmp;...
        eye(size(tmp,2))];
end

function B = formBspeedup(A_ie, A_ii, M_ee, M_ll, A_ln, A_el, lsolve)
    % Alternate B formulation with additional interpolation steps
    ll_ln = (lsolve(M_ll,A_ln));
    tmp = (lsolve(M_ee,A_el*ll_ln));
    B = [-(lsolve(A_ii,(A_ie*tmp)));...
        tmp;...
        -ll_ln;...
        eye(size(ll_ln,2))];
    return
end

function X = restrictOperator(CG, DG, Nf)
    %Returns a matrix representing the restriction operator
    %Should be nxf big where n is coarse nodes and f total fine nodes
    %Permute the partition to the new index space
    permuted_partition = DG.P*CG.partition;
    xind = zeros(1,Nf);
    yind = zeros(1,Nf);
    pos = 1;
    for coarse = 1:CG.cells.num
       %find indices of corresponding coarse nodes
       indices = find(permuted_partition == coarse);
       %insert the correct cells into the arrays of indices
       M = length(indices);
       yind(pos:(pos + M-1)) = indices;
       xind(pos:(pos + M-1)) = coarse;
       pos = pos + M;
    end
    X = sparse(xind, yind, 1) > 0;
    return
end

function Dmark = formD(Abar, CG, Nel)
    % Seven point stencil assumed for sparse allocation
    D = sparse([],[],[],Nel,Nel,7*Nel);
    index = 1;
    for coarse = 1:CG.cells.num
        tmp = sum(CG.partition == coarse);
        diagonal = index:(index - 1 + tmp);
        D(diagonal, diagonal) = Abar(diagonal, diagonal); %#ok fast enough
        index = index + tmp;
    end
    Dmark = D + diag(sum(Abar-D,1));
end

function [D U] = DU(Abar, CG, Nel)
    D = sparse([],[],[],Nel,Nel,7*Nel);
    U = sparse([],[],[],Nel,Nel,7*Nel);
    index = 1;
    for coarse = 1:CG.cells.num
        tmp = sum(CG.partition == coarse);
        diagonal = index:(index - 1 + tmp);
        D(diagonal, diagonal) = Abar(diagonal, diagonal); %#ok fast enough
        U(diagonal, diagonal(end)+1:end) = Abar(diagonal, diagonal(end)+1:end); %#ok
        index = index + tmp;
    end
end

% Flux for sides

function state = setFlux(G, p, T, fluid, state, opt)
    %COPIED FROM incompTPFA

   % Preliminaries
   cellNo = rldecode(1:G.cells.num, double(diff(G.cells.facePos)), 2).';
   cf     = G.cells.faces(:,1);
   nf     = G.faces.num;
   nc     = G.cells.num;
%    nw     = length(opt.wells);
%    n      = nc + nw;

   % Face transmissibility = harmonic average of half-transmissibilities
   [mob, omega, rho] = dynamic_quantities(state, fluid);
   totmob = sum(mob, 2);
   T      = T .* totmob(cellNo);
   ft     = 1 ./ accumarray(cf, 1./T, [nf, 1]);

   if(~isempty(opt.c2cTrans))
       [~,totmob_hyb] = computeHybridTrans(G,totmob(cellNo));
       ft_hyb 	= opt.c2cTrans.*totmob_hyb;
   end


   % Identify internal faces
   i  = all(G.faces.neighbors ~= 0, 2);

   % Boundary conditions and source terms.
   [ff, gg, hh, grav, dF, dC] = computePressureRHS(G,...
                                                   omega, ...
                                                   opt.bc,...
                                                   opt.src);

   sgn = 2*(G.faces.neighbors(cf, 1) == cellNo) - 1;
   j   = i(cf) | dF(cf);
   fg  = accumarray(cf(j), grav(j).*sgn(j), [nf, 1]);


   % Reconstruct face pressures and fluxes.
   fpress     =  ...
          accumarray(G.cells.faces(:,1), (p(cellNo)+grav).*T, [G.faces.num,1])./ ...
          accumarray(G.cells.faces(:,1), T, [G.faces.num,1]);


   % Neumann faces
   b         = any(G.faces.neighbors==0, 2);
   fpress(b) = fpress(b) - hh(b)./ft(b);


   % Contribution from gravity
   %fg         = accumarray(cf, grav.*sgn, [nf, 1]);
   %fpress(~i) = fpress(~i) + fg(~i);

   % Dirichlet faces
   fpress(dF) = dC;


   % Sign for boundary faces
   sgn  = 2*(G.faces.neighbors(~i,2)==0)-1;
   ni   = G.faces.neighbors(i,:);
   flux = -accumarray(find(i),  ft(i) .*(p(ni(:,2))-p(ni(:,1))-fg(i)), [nf, 1]);
   c    = sum(G.faces.neighbors(~i,:),2) ;
   fg  = accumarray(cf, grav, [nf, 1]);
   flux(~i) = -sgn.*ft(~i).*( fpress(~i) - p(c) - fg(~i) );
   %flux = -sgn.*ft((fpress(~i)-p(c)-grav));
   %state.pressure(1 : nc) = p(1 : nc);
   state.flux(:)          = flux;
   state.facePressure     = fpress;

   % Compute cell2cell flux
   if(~isempty(opt.c2cTrans))
       n_hyb = G.cells.neighbors;
       flux_hyb =-ft_hyb .*(p(n_hyb(:,2))-p(n_hyb(:,1)));
   else
       flux_hyb=[];
   end
   state.fluxc2c = flux_hyb;

   W = opt.wells;
   for k = 1 : numel(W),
      wc       = W(k).cells;
%       if strcmp(W(k).type,'bhp')
%           % Ok already handled
%           break;
%       end
      dp       = norm(gravity()) * W(k).dZ*sum(rho .* W(k).compi, 2);
      state.wellSol(k).flux     = W(k).WI.*totmob(wc).*(p(nc+k) + dp - p(wc));
      state.wellSol(k).pressure = p(nc + k);
   end
end

function [mob, omega, rho] = dynamic_quantities(state, fluid)
   [mu, rho] = fluid.properties(state);
   s         = fluid.saturation(state);
   kr        = fluid.relperm(s, state);

   mob    = bsxfun(@rdivide, kr, mu);
   totmob = sum(mob, 2);
   omega  = sum(bsxfun(@times, mob, rho), 2) ./ totmob;
end

%% Iterative helpers
function GinvU = Ginv(U, Ctimes, M_nn, Nf, R, X, A, B)

    CU = Ctimes(U);
    tmp = (speye(Nf) - R'*R)*U + R'*X*(U - A*CU);
    Ctmp = Ctimes(tmp);
    GinvU = (B*mldivide(M_nn,R*tmp) + Ctmp);
end



