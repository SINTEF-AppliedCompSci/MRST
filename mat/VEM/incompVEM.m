function state = incompVEM(state, G, S, fluid, varargin)
%Solve incompressible flow problem (fluxes/pressures) using first or second
%order virtual element method.
%
% SYNOPSIS:
%   state = incompVEM(state, G, S, fluid)
%   state = incompVEM(state, G, S, fluid, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function assembles and solves a set of linear equations defining
%   the pressure at the nodes (first order), faces, edges and cells (second
%   order) for the reservoir simulation problem defined by Darcy's law,
%   sources and boundary conditions.
%
% REQUIRED PARAMETERS:
%   state  - Reservoir and well solution structure either properly
%            initialized from function 'initState', or the results from a
%            previous call to function 'incompVEM' and, possibly, a
%            transport solver such as function 'explicitTransport'.
%
%   G, S   - Grid and (VEM) linear system data structures as defined by
%            function 'computeVirtualIP'.
%
%   fluid  - Fluid data structure as described by 'fluid_structure'.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   bc     - Boundary condition structure as defined by function 'addBC'.
%            This structure accounts for all external boundary conditions
%            to the reservoir flow.  May be empty (i.e., bc = []) which is
%            interpreted as all external no-flow (homogeneous Neumann)
%            conditions. Can also be a strucutre defined by the function
%            'addBCFunc', wich allows for function handle boundary
%            conditions for easy patch testing.
%
%   src    - Explicit source contributions as defined by function
%            'addSource'.  May be empty (i.e., src = []) which is
%            interpreted as a reservoir model without explicit sources.
%
%   facePressure -
%            Whether or not to calculate face pressures if a first order
%            method is used. Defalut value: facePressure = FALSE.
%
%   cellPressure -
%            Whether or not to calculate cell pressures if a first order
%            method is used. Defalut value: cellPressure = FALSE.
%
%   conservativeFlux - 
%           Wheter or not to apply postprocessing to make the calculated
%           flux field conservative in the local sense, necessary on order
%           to use incompVEM in transport solvers. Defalut value:
%           conservativeFlux = FALSE.
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
%            caller of function 'incompMimetic'.
%            Logical.  Default value: MatrixOutput = FALSE.
%
% RETURNS:
%   state - Update reservoir and well solution structure with new values
%           for the fields:
%              - nodePressure -- Pressure values for all nodes in the
%                                discretized resrvoir model, 'G'.
%              - edgePressure -- If G.griddim = 3, pressure values for all
%                                edges in the discretized resrvoir model,
%                                'G'.
%              - facePressure -- If G.griddim = 3, pressure values for all
%                                faces in the discretized resrvoir model,
%                                'G'.
%              - pressure     -- Pressure values for all cells in the
%                                discretised reservoir model, 'G'.
%              - flux         -- Flux across global interfaces
%                                corresponding to the rows of
%                                'G.faces.neighbors'.
%              - A        -- System matrix.  Only returned if specifically
%                            requested by setting option 'MatrixOutput'.
%
%
% NOTE:
%   If there are no external influences, i.e., if all of the structures
%   'W', 'bc', and 'src' are empty and there are no effects of gravity, and
%   no system right hand side has been supplied externally, then the input
%   state is returned unchanged and a warning is printed in the command
%   window.  This warning is printed with message ID
%
%           'incompMimetic:DrivingForce:Missing'
%
% SEE ALSO:
%   computeVirtualIP, addBC, addBCFunc addSource, initSimpleFluid
%   initState, solveIncompFlowMS.

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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

%   Written by Øystein Strengehagen Klemetsdal, SINTEF ICT/NTNU, 2016.

%   Merge input parameters
opt = struct('wells'           , []       , ...
             'bc'              , []       , ...
             'src'             , []       , ...
             'facePressure'    , false    , ...
             'cellPressure'    , false    , ...
             'conservativeFlux', false    , ...
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
                
    %   Assemble global matrix
    [A, P] = glob(state, G, S, fluid, N, tmob);
    
    %   Compute right-hand side
    rhs = computeRHS(G, S, P, N, opt);
    
    %   Add well contributions
    [A, rhs] = addWells(A, rhs, G, S, tmob, N, opt);
    
    
    
    %   Apply boundary condtitions.
    if ~isempty(opt.bc)
        [A, rhs] = imposeBC(A, rhs, G, S, opt.bc, N);
    end
    
    %   If no pressure bc is given, we follow the MRST convention of
    %   normalizing the pressure to zero in the first cell (or, for k = 1,
    %   the first node of the first cell.
    pressure_bc = ~isempty(opt.bc) && any(strcmpi('pressure', opt.bc.type));
    if ~pressure_bc && S.order == 2
%         if S.order == 1
%             n = G.cells.nodes(G.cells.nodePos(1):G.cells.nodePos(2)-1);
%             i = n(1);
%         else
            if G.griddim == 2  
                i = G.nodes.num + G.faces.num + 1;
            else
                i = G.nodes.num + G.edges.num + G.faces.num + 1;
            end
%         end
            A(i,i) = 2*A(i,i);
    end
    
end

%--------------------------------------------------------------------------

function [A, P] = glob(state, G, S, fluid, N, tmob)
%   Assemble global system.    

    %   Number of dofs per cell.
    NP = diff(G.cells.nodePos) ...
       + diff(G.cells.facePos)*polyDim(S.order-2, G.griddim-1) ...
       + polyDim(S.order-2, G.griddim);
    if G.griddim == 3
        NP = NP + diff(G.cells.edgePos)*polyDim(S.order-2, G.griddim-2);
    end
    
    %   Multiply by total mobility for each cell.
%     tmob = totmob(state, fluid);
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
            %   from teh def of the local VEM space.

            rhs = zeros(G.cells.num,1);
            rhs(src.cell) = src.rate;
            rhs = rldecode(rhs, diff(G.cells.nodePos), 1);
            PiNstar = squeezeBlockDiag(S.PiNstar, diff(G.cells.nodePos), ...
                 polyDim(S.order, G.griddim), sum(diff(G.cells.nodePos)))';
            rhs = rhs.*PiNstar(:,1);            
            rhs = P'*rhs;

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

function [A, rhs] = addWells(A, rhs, G, S, tmob, N, opt)

   W = opt.wells;
   
    
    if ~isempty(W),
          
        wc    = vertcat(W.cells);
        dofVec = G.nodes.num + G.edges.num*polyDim(S.order-2, G.griddim-2) ...
                             + G.faces.num*polyDim(S.order-2, G.griddim-1) ...
                             + wc;
        rhs(dofVec) = vertcat(W.val);
        v = ones(N,1);
%         v(dofVec) = tmob(wc).*vertcat(W.WI);
        I = spdiags(v, 0, N,N);
        A(dofVec, :) = I(dofVec,:);
        
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
        
        isNeu = find(strcmp(bc.type, 'flux'));
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

function p = facePressure(state, G, S, f)
     
    if G.griddim == 2
        
        n = G.faces.nodes(mcolon(G.faces.nodePos(f),G.faces.nodePos(f+1)-1));
        [ii, jj] = blockDiagIndex(ones(numel(f),1), 2*ones(numel(f),1));
        p = 1/2*sparse(ii, jj, 1)*state.nodePressure(n);

    else
        
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

    tm = totmob(state, fluid);
    
    ii    = [(1:G.cells.num)'; max(G.faces.neighbors, [], 2)];

    totmobMat = spdiags(tm(ii), 0, numel(ii),numel(ii));
    T      = S.T.T * totmobMat;
    
    bf = boundaryFaces(G);
    
    ii = [(1:G.cells.num)'; G.cells.num + bf];
    p = [state.pressure; state.facePressure(bf)];
    flux = cellFlux2faceFlux(G, T(:, ii)*p);
    
    neu = false(G.faces.num, 1);
    bf = boundaryFaces(G);
    neu(bf) = true;
    v = zeros(G.faces.num,1);
    if ~isempty(bc)
        neu(bc.face(strcmp(bc.type, 'pressure'))) = false;
        isNeu = strcmp(bc.type, 'flux');
        f = bc.face(isNeu);
        fSgn  = - 1 + 2*(G.faces.neighbors(f,2) ~= 0);
        v(f) = bc.value(isNeu).*fSgn;
    end
    flux(neu) = v(neu);
%     
%     flux(neu) = 0;
%     
%     if ~isempty(bc) && ~isfield(bc, 'func')
%         
%         isNeu = strcmp(bc.type, 'flux');
%         f = bc.face(isNeu);
%         fSgn  = - 1 + 2*(G.faces.neighbors(f,2) ~= 0);
%         flux(f) = bc.value(isNeu).*fSgn;
% 
%     end

end



function state = packSolution(state, p, G, S, fluid, opt, A)

[~, rho] = fluid.properties(state);

gvec = gravity();
if G.griddim == 3
    if S.order == 1
        pot = rho*gvec(3)*G.nodes.coords(:,3);
    else
        pot = rho*gvec(3)*[G.nodes.coords(:,3)   ; G.edges.centroids(:,3); ...
                           G.faces.centroids(:,3); G.cells.centroids(:,3)];
    end
    p = p-pot;
end

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
                       
if opt.facePressure
    state.facePressure = facePressure(state, G, S, 1:G.faces.num);
end

if (opt.cellPressure || ~isempty(S.T)) && isempty(state.pressure)
    state.pressure = cellPressure(state, G, S);
    if isempty(state.facePressure)
        state.facePressure(boundaryFaces(G)) = facePressure(state, G, S, boundaryFaces(G));
    end
end
    
if ~isempty(S.T)
    state.flux = computeFlux(state, G, S, fluid, opt.bc);
end

if opt.matrixOutput
    state.A   = A;
end

end

%--------------------------------------------------------------------------

function tm = totmob(state, fluid)

   [mu, ~] = fluid.properties(state);
   s       = fluid.saturation(state);
   kr      = fluid.relperm(s, state);
   
   mob    = bsxfun(@rdivide, kr, mu);
   tm = sum(mob, 2);

end

%--------------------------------------------------------------------------

function nk = polyDim(k, dim)

    if k < 0
        nk = 0;
    else
        nk = nchoosek(k+dim,k);
    end
end
