function state = incompVEM(state, G, S, fluid, varargin)
%Solve incompressible flow problem (fluxes/pressures) using a first- or
%second-order virtual element method.
%
% SYNOPSIS:
%   state = incompVEM(state, G, S, fluid)
%   state = incompVEM(state, G, S, fluid, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function assembles and solves a set of linear equations defining
%   the pressure at the nodes (first order), faces, edges and cells (second
%   order) for the reservoir simulation problem defined by Darcy's law,
%   sources and boundary conditions. The fluxes are reconstructed from the
%   pressure solution.
%
% REQUIRED PARAMETERS:
%   state  - Reservoir and well solution structure either properly
%            initialized from function 'initState', or as the results from
%            a previous call to function 'incompVEM' and, possibly, a
%            transport solver such as function 'explicitTransport'.
%
%   G, S   - Grid and (VEM) linear system data structures as defined by
%            function 'computeVirtualIP'.
%
%   fluid  - Fluid data structure as described by 'fluid_structure'.
%
% OPTIONAL PARAMETERS:
%   bc     - Boundary condition structure as defined by function 'addBC'.
%            This structure accounts for all external boundary conditions
%            to the reservoir flow.  May be empty (i.e., bc = []) which is
%            interpreted as all external no-flow (homogeneous Neumann)
%            conditions. Can also be a strucutre defined by the function
%            'addBCVEM', wich allows for function handle boundary
%            conditions for easy patch testing.
%
%   src    - Explicit source contributions as defined by function
%            'addSource'.  May be empty (i.e., src = []) which is
%            interpreted as a reservoir model without explicit sources.
%
%   facePressure -
%            Whether or not to calculate face pressures if a first-order
%            method is used. Defalut value: facePressure = FALSE.
%
%   cellPressure -
%            Whether or not to calculate cell pressures if a first-order
%            method is used. Defalut value: cellPressure = FALSE.
%
%   LinSolve -
%            Handle to linear system solver software to which the fully
%            assembled system of linear equations will be passed.  Assumed
%            to support the syntax
%
%                        x = LinSolve(A, b)
%
%            in order to solve a system Ax=b of linear equations.
%            Default value: LinSolve = @mldivide (backslash).
%
%   MatrixOutput -
%            Whether or not to return the final system matrix 'A' to the
%            caller of function 'incompVEM'.
%            Logical.  Default value: MatrixOutput = FALSE.
%
% RETURNS:
%   state - Update reservoir solution structure with new values for the
%           fields:
%              - nodePressure -- Pressure values for all nodes in the
%                                discretized resrvoir model, 'G'.
%              - edgePressure -- If G.griddim = 3 and method order = 2,
%                                pressure values for all edges in the
%                                discretized resrvoir model, 'G'.
%              - facePressure -- If method order = 2 or facePressure =
%                                true, pressure values for all faces in the
%                                discretized resrvoir model, 'G'.
%              - pressure     -- If method order = 2 or cellPressure =
%                                true, Pressure values for all cells in the
%                                discretised reservoir model, 'G'.
%              - flux         -- If calculateFlux = true, fluxes across
%                                global interfaces corresponding to the
%                                rows of 'G.faces.neighbors'.
%              - A            -- System matrix.  Only returned if
%                                specifically requested by setting option
%                                'MatrixOutput'.
%
% SEE ALSO:
%   computeVirtualIP, addBC, addBCVEM addSource, initSimpleFluid initState.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

%   Written by Øystein Strengehagen Klemetsdal, SINTEF/NTNU, 2016.

%%  Main function
%--------------------------------------------------------------------------

%   Merge input parameters
opt = struct('bc'              , []       , ...
             'src'             , []       , ...
             'facePressure'    , false    , ...
             'cellPressure'    , false    , ...
             'calculateFlux'   , true     , ...
             'linSolve'        , @mldivide, ...
             'matrixOutput'    , false         );
             
opt = merge_options(opt, varargin{:});

%   Check options
assert(G.griddim == 2 || G.griddim == 3, 'Physical dimension must be 2 or 3.');

%   Assemble and solve global system
[A, rhs] = assembleSystem(state, G, S, fluid, opt);
p        = opt.linSolve(A, rhs);

%   Pack solution struct
state = packSolution(state, p, G, S, fluid, opt, A);

end

%%  Helper functions
%--------------------------------------------------------------------------

function [A, rhs] = assembleSystem(state, G, S, fluid, opt)
%   Assembles the global system by summing the local matrices and applying
%   boundary conditions.

    if G.griddim == 2
        G.edges.num = 0;
    end

    %   Total number of dofs.
    N = G.nodes.num + G.edges.num*polyDim(S.order - 2, G.griddim - 2) ...
                    + G.faces.num*polyDim(S.order - 2, G.griddim - 1) ...
                    + G.cells.num*polyDim(S.order - 2, G.griddim    );
    
    %   calculate total mobilities-
    tmob = totmob(state, fluid);
                
    %   Assemble global matrix A. P is a map from local to global dofs.
    [A, P] = glob(G, S, N, tmob);
    
    %   Compute right-hand side
    rhs = computeRHS(G, S, P, N, opt);
    
    %   Apply boundary conditions.
    if ~isempty(opt.bc)
        [A, rhs] = imposeBC(A, rhs, G, S, opt.bc, N);
    end
    
    %   If no pressure bc is given, we follow the MRST convention of
    %   normalizing the pressure to zero in the first cell (or, for k = 1,
    %   the first node of the first cell.
    pressure_bc = ~isempty(opt.bc) && any(strcmpi('pressure', opt.bc.type));
    if ~pressure_bc
        if S.order == 1
            n = G.cells.nodes(G.cells.nodePos(1):G.cells.nodePos(2)-1);
            i = n(1);
        elseif G.griddim == 2  
            i = G.nodes.num + G.faces.num + 1;
        else
            i = G.nodes.num + G.edges.num + G.faces.num + 1;
        end
        A(i,i) = 2*A(i,i);
    end
    
end

%--------------------------------------------------------------------------

function [A, P] = glob(G, S, N, tmob)
%   Assemble global system.    

    %   Number of dofs per cell.
    NP = diff(G.cells.nodePos) ...
       + diff(G.cells.facePos)*polyDim(S.order-2, G.griddim-1) ...
       + polyDim(S.order-2, G.griddim);
    if G.griddim == 3
        NP = NP + diff(G.cells.edgePos)*polyDim(S.order-2, G.griddim-2);
    end
    
    %   Multiply by total mobility for each cell.
    tmob = spdiags(rldecode(tmob, NP, 1),0, sum(NP), sum(NP));
    
    %   Sum contribution from each cell.
    P = sparse(1:numel(S.dofVec), S.dofVec, 1, numel(S.dofVec), N);
    A = P'*(tmob*S.A)*P;

end

function rhs = computeRHS(G, S, P, N, opt)
%   Compute Right-hand side.

    src = opt.src;
    if ~isempty(src)
        if S.order == 1
            %   The first moment over the cell can be calculated exactly
            %   from the def of the local VEM space.

            rhs = zeros(G.cells.num,1);
            rhs(src.cell) = src.rate;
            rhs = rldecode(rhs, diff(G.cells.nodePos), 1);
            PiNstar = squeezeBlockDiag(S.PiNstar, diff(G.cells.nodePos), ...
                 polyDim(S.order, G.griddim), sum(diff(G.cells.nodePos)))';
            rhs = rhs.*PiNstar(:,1);            
            rhs = full(P'*rhs);

        else
            %   The first moment over the cell is now a degree of freedom.
            
            rhs = zeros(N,1);
            if G.griddim == 2
                ii = src.cell + G.nodes.num + G.faces.num;
            else
                ii = src.cell + G.nodes.num + G.edges.num + G.faces.num;
            end
              rhs(ii) = src.rate;
        end
        
    else
        rhs = zeros(N,1);
    end
    
end

%--------------------------------------------------------------------------

function [A, rhs] = imposeBC(A, rhs, G, S, bc, N)
%   Impose boundary conditions. Can either be function handles for each
%   cell (for pathc testing purposes only), or MRST grid structure.

    nfn = diff(G.faces.nodePos);
    
    if ~isfield(bc, 'func')
        %   Standard MRST BC struc.
        
        NF = diff(G.faces.nodePos);
        
        %   Neumann boundary faces.
        isNeu = strcmp(bc.type, 'flux');            
        f = bc.face(isNeu);
        v = bc.value(isNeu);
        n = G.faces.nodes(mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1));
        
        if G.griddim == 2
            %   VEM basis funciton are polynomials on each face, use
            %   centroid rule (k=1) or three-point (k=2) Gauss-Lobatto.
            
            P = sparse(n, 1:numel(n), 1, N, numel(n));
            vn = rldecode(v, NF(f), 1);
            if S.order == 1
                rhs = rhs + 1/2*P*vn;
            else
                rhs = rhs + 1/6*P*vn;
                rhs(f + G.nodes.num) = 2/3*v;
            end

            
        else            
            %   For k = 1, we use the projection onto the linear monomials
            %   and the centroid rule to evaluate (g, \phi^i)_{0,f}.
            %   For k = 2, they are given from the degrees of freedom.
            
            if S.order == 1

                PiNFstar = squeezeBlockDiag(S.PiNFstar, ...
                             NF(G.cells.faces(:,1)), polyDim(S.order, G.griddim-1), ...
                             sum(NF(G.cells.faces(:,1))));
                pos = [0; cumsum(NF)] + 1;
                v = rldecode(v, NF(f), 1).*PiNFstar(1,mcolon(pos(f), pos(f+1)-1))';

                NFf = NF(f);
                ii = rldecode((1:numel(f))', NFf, 1);
                dof = [0; cumsum(NFf(1:end-1))] + 1;
                iiN = mcolon(dof, dof + nfn(f) - 1);
                fDof(iiN) = n;
                v = sum(sparse(fDof, ii, v, N, numel(f)),2);            
                rhs = rhs + v;

            else
                
                rhs(f + G.nodes.num + G.edges.num) = ...
                                    rhs(f + G.nodes.num + G.edges.num) + v;

            end

        end
        
        %   Dirichlet boundary faces. We assue the pressure is constant on
        %   each face.
        isDir = strcmp(bc.type, 'pressure');
        f = bc.face(isDir);
        n = G.faces.nodes(mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1));
        v = bc.value(isDir);
        
        rhs(n) = rldecode(v, NF(f), 1);

        if S.order == 1
            dofVec = n;
        else
            
            if G.griddim == 2
                
                rhs(f + G.nodes.num) = v;
                dofVec = [n; f + G.nodes.num];
            
            else
                
                nfe = diff(G.faces.edgePos);
                e = G.faces.edges(mcolon(G.faces.edgePos(f), ...
                                         G.faces.edgePos(f+1)-1));
                rhs(e + G.nodes.num) = rldecode(v, nfe(f), 1);
                rhs(f + G.nodes.num + G.edges.num) = v;
                dofVec = [n; e + G.nodes.num; f + G.nodes.num + G.edges.num];
                
            end
            
        end

        I = speye(N);
        A(dofVec,:) = I(dofVec, :);

    else
        %   VEM-specific BCs, given as function handels.
        
        %   Neumann faces, only supported in 2D.
        isNeu = find(strcmp(bc.type, 'flux'));
        if nnz(isNeu) > 0
        
            assert(G.griddim == 2, ...
                   'Neumann conditions for addBCVEM only supported in 2D');
        
        end
        
        for i = 1:numel(isNeu)
            
            f = bc.face{isNeu(i)};
            g = bc.func{isNeu(i)};
            n = G.faces.nodes(mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1));
            
            if G.griddim == 2
                
                if S.order == 1
                
                    v = bsxfun(@times, bsxfun(@plus, ...
                                reshape(g(G.nodes.coords(n,:)), 2, [])'/6, ...
                                g(G.faces.centroids(f,:))/3), G.faces.areas(f));
                    v = reshape(v', [], 1);
                    rhs = rhs + sparse(n, 1:numel(n), 1, N, numel(n))*v;
                    
                else
                    
                    vn = bsxfun(@times, g(G.nodes.coords(n,:))/6, ...
                                rldecode(G.faces.areas(f), nfn(f), 1));
                    vn = reshape(vn', [], 1);
                    vf = g(G.faces.centroids(f,:))*2/3.*G.faces.areas(f);
                    
                    v = sparse(n, 1:numel(n), 1, N, numel(n))*vn;
                    v(f + G.nodes.num) = vf;
                    rhs = rhs + v;
                    
                end
                
            else
                
            end

        end
        
        isDir = find(strcmp(bc.type, 'pressure'));
                    
        I = speye(N);
        
        for i = 1:numel(isDir)
            
            f = bc.face{isDir(i)};
            g = bc.func{isDir(i)};
            n = G.faces.nodes(mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1));
            n = unique(n);
            
            if S.order == 1
            
                dofVec = n;
                rhs(dofVec) = g(G.nodes.coords(n,:));

            else
                
                if G.griddim == 2
                
                    dofVec = [n; f + G.nodes.num];
                    rhs(dofVec) = [g(G.nodes.coords(n,:)); g(G.faces.centroids(f,:))];
                    
                else
                    
                    e = G.faces.edges(mcolon(G.faces.edgePos(f), G.faces.edgePos(f+1)-1));
                    dofVec = [n; e + G.nodes.num; f + G.nodes.num + G.edges.num];
                    rhs(dofVec) = [g(G.nodes.coords(n,:)); ...
                                   g(G.edges.centroids(e,:)); ...
                                   polygonInt3D(G, f, g, 3)./G.faces.areas(f)];
                               
                end
                
            end
            
            A(dofVec,:) = I(dofVec,:);

        end

    end
end

%--------------------------------------------------------------------------

function state = packSolution(state, p, G, S, fluid, opt, A)
% Packs state solution struct, and adds face- cell pressures and flux if
% requested.

    if G.griddim == 2
        G.edges.num = 0;
    end

    state.nodePressure ...
               = full(p(1:G.nodes.num));
    state.edgePressure ...
               = full(p((1:G.edges.num*polyDim(S.order-2, G.griddim-2)) ...
                           + G.nodes.num));
    state.facePressure ...
               = full(p((1:G.faces.num*polyDim(S.order-2, G.griddim-1)) ...
                           + G.nodes.num ...
                           + G.edges.num*polyDim(S.order-2, G.griddim-2)));
    state.pressure     ...
               = full(p((1:G.cells.num*polyDim(S.order-2, G.griddim  )) ...
                           + G.nodes.num ...
                           + G.edges.num*polyDim(S.order-2, G.griddim-2)...
                           + G.faces.num*polyDim(S.order-2, G.griddim-1)));
    
    %   Calculate face pressures.
    if opt.facePressure
        state.facePressure = facePressure(state, G, S, 1:G.faces.num);
    end

    %   Calculate cell pressures. This must also be done it we use
    %   two-point or multipoint transmissibilities to reconstruct the flux.
    if (opt.cellPressure || ~isempty(S.T)) && isempty(state.pressure)
        state.pressure = cellPressure(state, G, S);
        if isempty(state.facePressure)
            state.facePressure(boundaryFaces(G)) ...
                             = facePressure(state, G, S, boundaryFaces(G));
        end
    end
    
    %   Calculate flux.
    if opt.calculateFlux
        state.flux = computeFlux(state, G, S, fluid, opt.bc);
    end

    %   return global matrix.
    if opt.matrixOutput
        state.A = A;
    end

end

function p = facePressure(state, G, S, f)
%   Calculate face pressures.
     
    if G.griddim == 2
        %   The pressure is linear on the faces, and we use linear
        %   interpolation to calculate the face pressures.
        
        n = G.faces.nodes(mcolon(G.faces.nodePos(f),G.faces.nodePos(f+1)-1));
        [ii, jj] = blockDiagIndex(ones(numel(f),1), 2*ones(numel(f),1));
        p = 1/2*sparse(ii, jj, 1)*state.nodePressure(n);

    else
        %   We know that the face pressure is a function in the 2D VEM
        %   space on each face, and we use this and the projection operator
        %   \Pi^\nabla to calculate the face pressures.
        %   NOTE: This calculates the pressure at ALL faces...
        
        cf = G.cells.faces(:,1);
        nfn = diff(G.faces.nodePos);
        c = squeezeBlockDiag(S.PiNFstar, nfn(cf),...
              polyDim(S.order, G.griddim-1), sum(cf));
        c = c(1,:)';

        ii = rldecode((1:numel(cf))', nfn(cf), 1);
        jj = 1:sum(nfn(cf));
        I = sparse(ii,jj,1);
                
        e = G.faces.edges(mcolon(G.faces.edgePos(cf), G.faces.edgePos(cf+1)-1));
        n = G.edges.nodes(mcolon(G.edges.nodePos(e),G.edges.nodePos(e+1)-1));
        n = reshape(n, 2, [])';
        n(G.faces.edgeSign(f) == -1,:) = n(G.faces.edgeSign(f) == -1, 2:-1:1);
        n = n(:,1);
        
        p = full(I*(c.*state.nodePressure(n)));
        
        ii = cf;
        jj = 1:numel(cf);
        I = sparse(ii, jj, 1);
        p = I*p./sum(I,2);        
        p = p(f);
        
    end
end

function p = cellPressure(state, G, S)
%   Calculates the cell pressures using the node pressures and the
%   projection operator \Pi^\nabla.
    
    ncn = diff(G.cells.nodePos);
    c = squeezeBlockDiag(S.PiNstar, ncn, polyDim(S.order, G.griddim), sum(ncn));
    c = c(1,:)';
    ii = rldecode((1:G.cells.num)', ncn, 1);
    jj = 1:numel(G.cells.nodes);
    I = sparse(ii,jj,1);
    p = full(I*(c.*state.nodePressure(G.cells.nodes)));
    
end

%--------------------------------------------------------------------------

function flux = computeFlux(state, G, S, fluid, bc)
%   Reconstruct fluxes from pressure field. Can be done using TPFA of MPFA
%   scheme.

    tm = totmob(state, fluid);
    
    if isempty(S.T)
        
        f    = G.cells.faces(:,1);
        fSgn = (-ones(numel(f),1)).^(G.faces.neighbors(f,1) ...
               ~= rldecode((1:G.cells.num)', diff(G.cells.facePos), 1)); 
        if size(f,1) == 1; f = f'; end
        
        k = S.order;
        d = G.griddim;
        ncn = diff(G.cells.nodePos);
        ncf = diff(G.cells.facePos);
        if G.griddim == 2
            nce = zeros(G.cells.num,1);
        else
            nce = diff(G.cells.edgePos);
        end
        
        %   Dofs per cell
        NP = ncn + nce.*polyDim(k-2,d-2) ...
                 + ncf*polyDim(k-2, d-1) ...
                 + polyDim(k-2, d);
             
        %   Map from global to local node dofs.
        vec = [1; cumsum(NP(1:end-1)) + 1];
        iiN = mcolon(vec, vec + ncn-1);
        
        %   Initialize pressure vector and r = dof([x,y,z]).
        p = zeros(sum(NP),1);
        r = zeros(G.griddim, sum(NP));
        
        if G.griddim == 2
            
            %   Map from global to local face and cell dofs.
            iiF = mcolon(vec + ncn, vec + ncn + ncf*polyDim(k-2, d-2) -1);
            iiP = mcolon(vec + ncn + ncf*polyDim(k-2, d-2), ...
                         vec + ncn + ncf*polyDim(k-2, d-2) + polyDim(k-2, d-1) -1);
            
            %   Nodes for each face of each cell.
            n   = G.faces.nodes(mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1));
            if size(n,1) == 1; n = n'; end
            n   = reshape(n,2,[])';
            n(fSgn == -1,:) = n(fSgn == -1,2:-1:1);
            n   = n(:,1);
            
            %   Shuffle presure and polynomial dofs.
            if k == 1
                p([iiN, iiF, iiP]) = state.nodePressure(n);
                r(:, [iiN, iiF, iiP]) = G.nodes.coords(n,:)';
            else
                p([iiN, iiF, iiP]) = [state.nodePressure(n); ...
                                    state.facePressure(f); state.pressure];
                r(:, [iiN, iiF, iiP]) = [G.nodes.coords(n,:); ...
                               G.faces.centroids(f,:); G.cells.centroids]';
            end
            
        else
            
            %   Map from global to local edge, face and cell dofs.
            iiE = mcolon(vec + ncn, vec + ncn + nce*polyDim(k-2, 1) - 1);
            iiF = mcolon(vec + ncn + nce*polyDim(k-2, 1), ...
                         vec + ncn + nce*polyDim(k-2, 1) ...
                                   + ncf*polyDim(k-2, 2) - 1);
            iiP = mcolon(vec + ncn + nce*polyDim(k-2, 1) ...
                 + ncf*polyDim(k-2, 2), ...
                         vec + ncn + nce*polyDim(k-2, 1) ...
                                   + ncf*polyDim(k-2, 2) + polyDim(k-2,3)-1);
            
            %   Suffle pressure and polyomial dofs.
            if k == 1
                p([iiN, iiE, iiF, iiP]) = state.nodePressure(G.cells.nodes);
                r(:, [iiN, iiE, iiF, iiP]) = G.nodes.coords(G.cells.nodes,:)';
            else
                p([iiN, iiE, iiF, iiP]) ...
                    = [state.nodePressure(G.cells.nodes); ...
                       state.edgePressure(G.cells.edges); ...
                       state.facePressure(f); ...
                       state.pressure];
                r(:, [iiN, iiE, iiF, iiP]) ...
                                  = [G.nodes.coords(G.cells.nodes,:); ...
                                     G.edges.centroids(G.cells.edges,:); ...
                                     G.faces.centroids(f,:); ...
                                     G.cells.centroids]';
            end
        end
        
        %   Put in block matrix form
        [ii,jj] = blockDiagIndex(ones(G.cells.num,1), NP);
        pm = sparse(ii,jj,p');
        [jj,ii] = blockDiagIndex(G.griddim*ones(G.cells.num,1), NP);
        rm = sparse(ii, jj, r(:));
        
        %   Divide local matrices by total mobility.
        totmobMat = spdiags(rldecode(tm, NP, 1), 0, sum(NP), sum(NP));
        A = totmobMat*S.A;
        
        %   Calculate K \nabla p for each cell.
        Kgradp = bsxfun(@rdivide, squeezeBlockDiag((pm*A*rm), ...
            ones(G.cells.num,1), G.cells.num, G.griddim), G.cells.volumes);
        
        %   Calculate flux for each face of each cell.
        flux = sum(rldecode(-Kgradp, diff(G.cells.facePos),1)...
                                                 .*G.faces.normals(f,:),2);
        
        %   Compute face fluxes by arithmetic mean.
        P = sparse(f, 1:numel(f),1);
        P = bsxfun(@rdivide, P, sum(P,2));
        flux = P*flux;
        
    %   MPFA.
    elseif strcmp(S.transType, 'mpfa')
    
        ii    = [(1:G.cells.num)'; max(G.faces.neighbors, [], 2)];

        totmobMat = spdiags(tm(ii), 0, numel(ii),numel(ii));
        T      = S.T.T * totmobMat;

        bf = boundaryFaces(G);
        ii = [(1:G.cells.num)'; G.cells.num + bf];
        p = [state.pressure; state.facePressure(bf)];
        flux = cellFlux2faceFlux(G, T(:, ii)*p);

    %   TPFA.
    elseif strcmp(S.transType, 'tpfa')
        
        [neighborship, ~] = getNeighbourship(G, 'Topological', true);
        [cellNo, cf, ~] = getCellNoFaces(G);
        nif    = size(neighborship, 1);
              
        ii  = all(neighborship ~= 0, 2);
        ni   = neighborship(ii,:);
        
        T  = S.T .* tm(cellNo);
        ft = 1 ./ accumarray(cf, 1 ./ T, [nif, 1]);
        
        p = state.pressure;
        flux = -accumarray(find(ii),  ft(ii) .*(p(ni(:,2))-p(ni(:,1))), [nif, 1]);
        sgn  = 2*(G.faces.neighbors(~ii,2)==0)-1;
        c    = sum(G.faces.neighbors(~ii,:),2) ;
        flux(~ii) = -sgn.*ft(~ii).*( state.facePressure(~ii) - p(c));  
        
    end
    
    %   Replace flux on neumann faces by prescribed value. Not supported
    %   for addBCVEM.
    
    neu = false(G.faces.num, 1);
    bf = boundaryFaces(G);

    neu(bf) = true;
    v = zeros(G.faces.num,1);
    if ~isempty(bc) && ~isfield(bc, 'func')
        neu(bc.face(strcmp(bc.type, 'pressure'))) = false;
        isNeu = strcmp(bc.type, 'flux');
        f = bc.face(isNeu);
        fSgn  = - 1 + 2*(G.faces.neighbors(f,2) ~= 0);
        v(f) = bc.value(isNeu).*fSgn;
    end
    flux(neu) = v(neu);

end
        
%--------------------------------------------------------------------------

function tm = totmob(state, fluid)
%   Calculate total mobility.

    [mu, ~] = fluid.properties(state);
    s       = fluid.saturation(state);
    kr      = fluid.relperm(s, state);

    mob = bsxfun(@rdivide, kr, mu);
    tm  = sum(mob, 2);

end
