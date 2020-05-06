function state = incompMPFA_DFM(state, g, T, fluid, varargin)
%Solve incompressible flow problem (fluxes/pressures) using MPFA-O method.
%
% SYNOPSIS:
%   state = incompMPFA_DFM(state, G, T, fluid)
%   state = incompMPFA_DFM(state, G, T, fluid, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function assembles and solves a (block) system of linear equations
%   defining interface fluxes and cell pressures at the next time step in a
%   sequential splitting scheme for the reservoir simulation problem
%   defined by Darcy's law and a given set of external influences (wells,
%   sources, and boundary conditions).
%
%   This function uses a multi-point flux approximation (MPFA) method with
%   minimal memory consumption within the constraints of operating on a
%   fully unstructured polyhedral grid structure.
%
%   The function is modified from the original incompMPFA (legacy version) to account for
%   non-neighbor cell-to-cell connections for hybrid cells.
%
% REQUIRED PARAMETERS:
%   state  - Reservoir and well solution structure either properly
%            initialized from functions 'initResSol' and 'initWellSol'
%            respectively, or the results from a previous call to function
%            'incompMPFA' (legacy version) and, possibly, a transport solver such as
%            function 'implicitTransport'.
%
%   G, T   - Grid and half-transmissibilities as computed by the function
%            'computeTrans'.
%
%   fluid  - Fluid object as defined by function 'initSimpleFluid'.
%
% OPTIONAL PARAMETERS:
%   wells  - Well structure as defined by functions 'addWell' and
%            'assembleWellSystem'.  May be empty (i.e., W = struct([]))
%            which is interpreted as a model without any wells.
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
%   LinSolve     - Handle to linear system solver software to which the
%                  fully assembled system of linear equations will be
%                  passed.  Assumed to support the syntax
%
%                        x = LinSolve(A, b)
%
%                  in order to solve a system Ax=b of linear equations.
%                  Default value: LinSolve = @mldivide (backslash).
%
%   MatrixOutput - Whether or not to return the final system matrix 'A' to
%                  the caller of function 'incompTPFA'.
%                  Logical.  Default value: MatrixOutput = FALSE.
%
%   Verbose      - Whether or not to time portions of and emit informational
%                  messages throughout the computational process.
%                  Logical.  Default value dependent on global verbose
%                  setting in function 'mrstVerbose'.
%
%   cellConnection -If true, cell connections spesified by T.cellTrans are added
%                   The connections must be spesified in G.cells.neighbors.
%                   A vector flux2 is returned contining the extra fluxes.
%
%   upwind         -If true, the upwind mobility is used, els the harmonic
%                   averange of the mobility is used
%
% RETURNS:
%   xr - Reservoir solution structure with new values for the fields:
%          - pressure     -- Pressure values for all cells in the
%                            discretised reservoir model, 'G'.
%          - boundaryPressure --
%                            Pressure values for all boundary interfaces in
%                            the discretised reservoir model, 'G'.
%          - flux         -- Flux across global interfaces corresponding to
%                            the rows of 'G.faces.neighbors'.
%          - A            -- System matrix.  Only returned if specifically
%                            requested by setting option 'MatrixOutput'.
%
%   xw - Well solution structure array, one element for each well in the
%        model, with new values for the fields:
%           - flux     -- Perforation fluxes through all perforations for
%                         corresponding well.  The fluxes are interpreted
%                         as injection fluxes, meaning positive values
%                         correspond to injection into reservoir while
%                         negative values mean production/extraction out of
%                         reservoir.
%           - pressure -- Well pressure.
%
% NOTE:
%   If there are no external influences, i.e., if all of the structures
%   'W', 'bc', and 'src' are empty and there are no effects of gravity,
%   then the input values 'xr' and 'xw' are returned unchanged and a
%   warning is printed in the command window. This warning is printed with
%   message ID
%
%           'incompMPFA:DrivingForce:Missing'
%
% EXAMPLE:
%    G   = computeGeometry(cartGrid([3,3,5]));
%    f   = initSingleFluid();
%    rock.perm = rand(G.cells.num, 1)*darcy()/100;
%    bc  = pside([], G, 'LEFT', 2);
%    src = addSource([], 1, 1);
%    W   = verticalWell([], G, rock, 1, G.cartDims(2), ...
%                       (1:G.cartDims(3)), 'Type', 'rate', 'Val', 1/day(), ...
%                       'InnerProduct', 'ip_tpf');
%    W   = verticalWell(W, G, rock, G.cartDims(1),   G.cartDims(2), ...
%                       (1:G.cartDims(3)), 'Type', 'bhp', ...
%                       'Val',  1*barsa(), 'InnerProduct', 'ip_tpf');
%    T   = computeMultiPointTrans(G, rock);
%    xr  = initResSol (G, 10);
%    xw  = initWellSol(G, 10);
%    [xr,xw] = incompMPFA(xr, xw, G, T, f, 'bc',bc,'src',src,'wells',W,...
%                         'MatrixOutput',true);
%
%    plotCellData(G, xr.cellPressure);
%
% SEE ALSO:
%   `computeMultiPointTransLegacy`, `addBC`, `addSource`, `addWell`, `initSingleFluid`,
%   `initResSol`, `initWellSol`, `mrstVerbose`.

%{
Copyright 2009, 2010, 2011 SINTEF ICT, Applied Mathematics.

Portions Copyright 2011-2012 University of Bergen.

This file is part of DFM module of The MATLAB Reservoir Simulation Toolbox
(MRST).

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


% Written by Jostein R. Natvig, SINTEF ICT, 2009.

   opt = struct('bc', [], 'src', [], 'wells', [], ...
                'LinSolve',     @mldivide,        ...
                'MatrixOutput', false,            ...
                'Verbose',mrstVerbose,'cellConnections',false,'upwind',false);
   opt = merge_options(opt, varargin{:});

   g_vec   = gravity();
   no_grav = ~(norm(g_vec) > 0); %(1 : size(g.nodes.coords,2))) > 0);
   if all([isempty(opt.bc)   , ...
           isempty(opt.src)  , ...
           isempty(opt.wells), no_grav]),
      warning(id('DrivingForce:Missing'),                       ...
              ['No external driving forces present in model--', ...
               'state remains unchanged.\n']);
   end

   cellNo = rldecode(1:g.cells.num, diff(g.cells.facePos), 2) .';
   cf     = g.cells.faces(:,1);
   nf     = g.faces.num;
   nc     = g.cells.num;
   nw     = length(opt.wells);
   n      = nc + nw;

   [totmob, omega, rho] = dynamic_quantities(state, fluid);

   % Needed after introduction of gravity
   Tg = T.Tg;
   T  = T.T;
   C = sparse(1:size(g.cells.faces, 1), cellNo, 1);

   % Computing the total transmissibility by multipling the geometric
   % transmissibility with the mobility.
   % If spesified the upwind mobility is used, els a harmonic average is used
   % for the mobilities.
   if(opt.upwind)

       % find the upwind index
       a = state.flux<0;
       a = a(cf);
       N = g.faces.neighbors(cf,:);
       intern = all(N~=0,2);
       upwind = zeros(length(N),1);
       upwind(~a & intern) = N(~a & intern,1);
       upwind(a & intern) = N(a & intern,2);
       upwind(~intern) = sum(N(~intern,:),2);

       % store the upwind mobility.
       totmob_upwind = totmob(upwind);

       % find the upwind mobility for the cell2cell connections.
       if(opt.cellConnections)

           % a map between an intersection its neigbor faces.
           cFace=[rldecode(1:length(g.hybridNeighbors.n),diff(g.hybridNeighbors.facePos),2)',g.faces.neighbors(g.hybridNeighbors.faces,1)];

           % the hybrid flux
           hflux=(state.flux(g.hybridNeighbors.faces));

           % the hybrid upwind index
           upwind2=hflux>0;

           % use a weighted upwind mobility:
           %
           %    mob_0 = sum ( (mob_j q_j)/sum(q_k))
           %
           % for all upwind cell j i.e q_j>0;
           %

           sumFlux=accumarray(cFace(upwind2>0,1),hflux(upwind2>0),[length(g.hybridNeighbors.n),1]);
           mob0=accumarray(cFace(upwind2,1),totmob(cFace(upwind2,2)).*hflux(upwind2),[length(g.hybridNeighbors.n),1])./sumFlux;

           % map from intersection to face
           mob0=mob0(cFace(:,1));

           % map from face to cellFace
           [~,cellFaces_hybrid]=ismember(g.hybridNeighbors.faces,cf);

           % the mobility is updated with the new upwind mobilities for the cell2cell connections.
           totmob_upwind(cellFaces_hybrid(~isnan(mob0) & ~upwind2))=mob0(~isnan(mob0) & ~upwind2);
       end

       % total transmissibility
       totmob_mat = spdiags(totmob_upwind , 0, numel(intern),numel(intern));
       T      =  totmob_mat*T;

   else
       % the original code, what happens?
       ind    = [(1:g.cells.num)'; max(g.faces.neighbors, [], 2)];
       T=T*spdiags(totmob(ind),0,numel(ind),numel(ind));

       % is used for the gravity
       totmob_mat = spdiags(rldecode(totmob, diff(g.cells.facePos)), 0, ...
                     size(g.cells.faces,1), size(g.cells.faces,1));
   end


   Tg     = totmob_mat*Tg;

   % identify internal faces
   i  = all(g.faces.neighbors ~= 0, 2);

   % Boundary conditions and source terms.
   % Note: Function 'computeRHS' is a modified version of function
   % 'computePressureRHS' which was originally written to support the
   % hybrid mimetic method.
   [ff, gg, hh, grav, dF, dC] = computePressureRHS(g, omega, ...
                                                   opt.bc, opt.src);
   b  = any(g.faces.neighbors==0, 2);
   if(opt.cellConnections)
       b(g.hybridNeighbors.faces)=false;
   end
   I1 = [(1:g.cells.num)'; g.cells.num + find(b)];
   D  = sparse(1:size(g.cells.faces,1), double(g.cells.faces(:,1)), 1);
   A  = [C, -D(:,b)]' * T(:,I1);

   % Gravity contribution for each face
   %cf  = g.cells.faces(:,1);
   %j   = i(cf) | dF(cf);
   %s   = 2*(g.faces.neighbors(cf, 1) == cellNo) - 1;
   fg  = [C, -D(:,b)]' * (Tg * grav);
   %fg  = accumarray(cf(j), grav(j).*s(j), [g.faces.num, 1]);

   rhs = [gg; -hh(b)];

   %%
   %if(opt.cellConnections)
   % A=A+sparse([N2(:,1),N2(:,1),N2(:,2),N2(:,2)],[N2(:,1),N2(:,2),N2(:,1),N2(:,2)],[-T2,T2,T2,-T2],length(I1),length(I1));
   %end

   %% Eliminate all but the cellpressure
   %BB=A(nc+1:end,nc+1:end);
   %AA=A(1:nc,1:nc);
   %DD=A(nc+1:end,1:nc);
   %DU=A(1:nc,nc+1:end);

   %A=AA-DU*inv(BB)*DD;
   %rhs=rhs(1:nc)+DU*inv(BB)*rhs(nc+1:end);

   %B=A(nc+1:end,
   %A=inv(A(nc+1:end)A(1:nc,1:nc)

   %% Dirichlet condition
   % If there are Dirichlet conditions, move contribution to rhs.  Replace
   % equations for the unknowns by speye(*)*x(dF) = dC.
   %% add gravity

   factor = A(1,1);
   assert (factor > 0)
   if any(dF),
      ind        = [false(g.cells.num, 1) ; dF(b)];
      rhs        = rhs - A(:,ind)*dC;
      rhs(ind)   = factor*dC;
      A(ind,:)   = 0;
      A(:,ind)   = 0;
      A(ind,ind) = factor * speye(sum(ind));
   end

   nnp=length(rhs);
   rhs=rhs-fg;
   rhs=[rhs;zeros(nw, 1)];

   %%%%%%%%%%%%%%%%%%%
   % add well equations
   C    = cell (nnp, 1);
   D    = zeros(nnp, 1);
   W    = opt.wells;
   d  = zeros(g.cells.num, 1);
   for k = 1 : nw,
      wc       = W(k).cells;
      nwc      = numel(wc);
      w        = k + nnp;

      wi       = W(k).WI .* totmob(wc);

      dp       = norm(gravity()) * W(k).dZ*sum(rho .* W(k).compi, 2);
      d   (wc) = d   (wc) + wi;

      if     strcmpi(W(k).type, 'bhp'),
         ww=max(wi);
         %ww=1.0;
         rhs (w)  = rhs (w)  + ww*W(k).val;
         rhs (wc) = rhs (wc) + wi.*(W(k).val + dp);
         C{k}     = -sparse(1, nnp);
         D(k)     = ww;

      elseif strcmpi(W(k).type, 'rate'),
         rhs (w)  = rhs (w)  + W(k).val;
         rhs (wc) = rhs (wc) + wi.*dp;

         C{k}     =-sparse(ones(nwc, 1), wc, wi, 1, nnp);
         D(k)     = sum(wi);

         rhs (w)  = rhs (w) - wi.'*dp;

      else
         error('Unsupported well type.');
      end
   end

   C = vertcat(C{:});
   D = spdiags(D, 0, nw, nw);
   A = [A, C'; C D];
   A = A+sparse(1:nc,1:nc,d,size(A,1),size(A,2));

   %if norm(gravity()) > 0,
%           rhs = rhs + T(:,I1)'*fg(cf);
%   end
   if ~any(dF) && (isempty(W) || ~any(strcmpi({W.type }, 'bhp'))),
      A(1) = A(1)*2;
   end
   ticif(opt.Verbose);
   p = opt.LinSolve(A, rhs);

   tocif(opt.Verbose);

   %% ---------------------------------------------------------------------
   dispif(opt.Verbose, 'Computing fluxes, face pressures etc...\t\t');
   ticif (opt.Verbose);

   % Reconstruct face pressures and fluxes.
   %fpress     =  ...
%          accumarray(g.cells.faces(:,1), (p(cellNo)+grav).*T, [g.faces.num,1])./ ...
%          accumarray(g.cells.faces(:,1), T, [G.faces.num,1]);


   % Neumann faces
   b         = any(g.faces.neighbors==0, 2);
   %fpress(b) = fpress(b) - hh(b)./ft(b);


   % Contribution from gravity
   %fg         = accumarray(cf, grav.*sgn, [nf, 1]);
   %fpress(~i) = fpress(~i) + fg(~i);

   % Dirichlet faces
   %fpress(dF) = dC;


   % Sign for boundary faces
   %sgn  = 2*(G.faces.neighbors(~i,2)==0)-1;
   %ni   = G.faces.neighbors(i,:);
   %flux = -accumarray(find(i),  ft(i) .*(p(ni(:,2))-p(ni(:,1))-fg(i)), [nf, 1]);
   %c    = sum(G.faces.neighbors(~i,:),2) ;
   %fg  = accumarray(cf, grav, [nf, 1]);
   %flux(~i) = -sgn.*ft(~i).*( fpress(~i) - p(c) - fg(~i) );
   %flux = -sgn.*ft((fpress(~i)-p(c)-grav));

   state.pressure(1 : nc) = p(1 : nc);
   state.flux(:)          = cellFlux2faceFlux(g, T(:, I1) * p(1 : nnp)+Tg * grav);

   % The cell2cell flux is weighed as following
   %
   %  q_ij = q_i0 q_j0 / sum_j q_j0
   %
   % For motivation and details, see Sandve, Berre, Nordbotten, JCP 2012
   if(opt.cellConnections)

       % the hybrid flux
       flux_h = state.flux(g.hybridNeighbors.faces);

       % the hybrid faces
       faces_h = rldecode(1:length(g.hybridNeighbors.n),diff(g.hybridNeighbors.facePos),2)';

       % number of hybrid segments meeting at an intersection
       n_f = rldecode(1:length(g.hybridNeighbors.n),g.hybridNeighbors.n,2)';

       % the sum of the incomming fluxes
       sumFlux = accumarray(faces_h(flux_h>0),flux_h(flux_h>0));

       % the weighted flux
       flux_ij = prod(flux_h(g.hybridNeighbors.neighbors),2)./sumFlux(n_f);

       % only the negative fluxes are stored
       flux_ij(isnan(flux_ij) | flux_ij>0) = 0;

       % the sign is changed such that the flux goes from neighor 1 to 2
       sgn = 2*(flux_h(g.hybridNeighbors.neighbors(:,1))<0)-1;

       % the cell2cell flux is stored with correct sign.
       state.fluxc2c = flux_ij.*sgn;
   else
       % if there is no cell2cell connections a empty array is returned.
       state.fluxc2c = [];
   end

   state.boundaryPressure = p(nc + 1 : nnp);

   for k = 1 : nw,
      wc       = W(k).cells;
      dp       = norm(gravity()) * W(k).dZ*sum(rho.*W(k).compi, 2);
      state.wellSol(k).flux = W(k).WI.*totmob(wc).*(p(nnp+k) + dp - p(wc));
      state.wellSol(k).pressure = p(nnp + k);
   end

   if opt.MatrixOutput,
      state.A   = A;
      state.rhs = rhs;
   end



   tocif(opt.Verbose);
end

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function s = id(s)
   s = ['incompMPFA:', s];
end

%--------------------------------------------------------------------------

function [totmob, omega, rho] = dynamic_quantities(state, fluid)
   [mu, rho] = fluid.properties(state);
   s         = fluid.saturation(state);
   kr        = fluid.relperm(s, state);

   mob    = bsxfun(@rdivide, kr, mu);
   totmob = sum(mob, 2);
   omega  = sum(bsxfun(@times, mob, rho), 2) ./ totmob;
end
