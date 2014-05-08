function state = incompTPFAVE_coupled(state, G, T, fluid, varargin)
%Solve incompressible flow problem (fluxes/pressures) using TPFA method.
%
% SYNOPSIS:
%   state = incompTPFA(state, G, T, fluid)
%   state = incompTPFA(state, G, T, fluid, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function assembles and solves a (block) system of linear equations
%   defining interface fluxes and cell pressures at the next time step in a
%   sequential splitting scheme for the reservoir simulation problem
%   defined by Darcy's law and a given set of external influences (wells,
%   sources, and boundary conditions).
%
%   This function uses a two-point flux approximation (TPFA) method with
%   minimal memory consumption within the constraints of operating on a
%   fully unstructured polyhedral grid structure.
%
% REQUIRED PARAMETERS:
%   state  - Reservoir and well solution structure either properly
%            initialized from functions 'initResSol' and 'initWellSol'
%            respectively, or the results from a previous call to function
%            'incompTPFA' and, possibly, a transport solver such as
%            function 'implicitTransport'.
%
%   G, T   - Grid and half-transmissibilities as computed by the function
%            'computeTrans'.
%
%   fluid  - Fluid object as defined by function 'initSimpleFluid'.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
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
%   verbose - Enable output.  Default value dependent upon global verbose
%             settings of function 'mrstVerbose'.
%
%   condition_number - Display estimated condition number of linear system.
%
% RETURNS:
%   state - Update reservoir and well solution structure with new values
%           for the fields:
%              - pressure -- Pressure values for all cells in the
%                            discretised reservoir model, 'G'.
%              - facePressure --
%                            Pressure values for all interfaces in the
%                            discretised reservoir model, 'G'.
%              - flux     -- Flux across global interfaces corresponding to
%                            the rows of 'G.faces.neighbors'.
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
%
% NOTE:
%   If there are no external influences, i.e., if all of the structures
%   'W', 'bc', and 'src' are empty and there are no effects of gravity,
%   then the input values 'xr' and 'xw' are returned unchanged and a
%   warning is printed in the command window. This warning is printed with
%   message ID
%
%           'incompTPFA:DrivingForce:Missing'
%
% EXAMPLE:
%    G   = computeGeometry(cartGrid([3, 3, 5]));
%
%    f   = initSingleFluid('mu' , 1*centi*poise, ...
%                          'rho', 1000*kilogram/meter^3);
%    rock.perm = rand(G.cells.num, 1)*darcy()/100;
%
%    bc  = pside([], G, 'LEFT', 2*barsa);
%    src = addSource([], 1, 1);
%    W   = verticalWell([], G, rock, 1, G.cartDims(2), []   , ...
%                       'Type', 'rate', 'Val', 1*meter^3/day, ...
%                       'InnerProduct', 'ip_tpf');
%    W   = verticalWell(W, G, rock, G.cartDims(1), G.cartDims(2), [], ...
%                       'Type', 'bhp', 'Val', 1*barsa, ...
%                       'InnerProduct', 'ip_tpf');
%
%    T   = computeTrans(G, rock);
%
%    state         = initResSol (G, 10);
%    state.wellSol = initWellSol(G, 10);
%
%    state = incompTPFA(state, G, T, f, 'bc', bc, 'src', src, ...
%                       'wells', W, 'MatrixOutput', true);
%
%    plotCellData(G, state.pressure)
%
% SEE ALSO:
%   computeTrans, addBC, addSource, addWell, initSingleFluid, initResSol,
%   initWellSol.

%{
#COPYRIGHT#
%}

% $Date: 2012-10-05 10:15:39 +0200 (Fri, 05 Oct 2012) $
% $Revision: 10006 $

   opt = struct('bc', [], 'src', [], 'wells', [], ...
                'LinSolve',     @mldivide,        ...
                'MatrixOutput', false, ...
                'Verbose',      mrstVerbose,...
                'condition_number',false,...
                'pc_form','wetting', ...
                'region3D', []);

   opt = merge_options(opt, varargin{:});

   g_vec   = gravity();
   if(any(strcmp(G.type, 'topSurfaceGrid')))
      no_grav=~(norm(gravity)>0);
   else
      no_grav = ~(norm(g_vec(1 : size(G.nodes.coords,2))) > 0);
   end
   if all([isempty(opt.bc)   , ...
           isempty(opt.src)  , ...
           isempty(opt.wells), no_grav]),
      warning(msgid('DrivingForce:Missing'),                   ...
             ['No external driving forces present in model--', ...
              'state remains unchanged.\n']);
   end

   %% ---------------------------------------------------------------------
   dispif(opt.Verbose, 'Setting up linear system...\t\t\t');
   ticif (opt.Verbose);

   bndFaces2D3D = G.facesBnd.index;
   % Preliminaries
   cellNo = rldecode(1:G.cells.num, double(diff(G.cells.facePos)), 2).';
   cf     = G.cells.faces(:,1);
   nf     = G.faces.num;
   nc     = G.cells.num;
   nw     = length(opt.wells);
   n      = nc + nw;
   
   % compute mobilities and gravity correction for cells on 2D bnd
   
   
   % Face transmissibility = harmonic average of half-transmissibilities
   % NB: omega is not correct!!
   [gravCorrVE, mob, mobCF, omegaCF, rho] = dynamic_quantities(state, fluid, G, cellNo, opt);   
   
   totmobCF = sum(mobCF, 2);
   T      = T .* totmobCF; 
   ft     = 1 ./ accumarray(cf, 1./T, [nf, 1]);  
   
   % Identify internal faces
   i  = all(G.faces.neighbors ~= 0, 2);

   % Boundary conditions and source terms.
   
   [ff, gg, hh, grav, dF, dC] = computePressureRHSVE_coupled(G,...
                                                   omegaCF, ...
                                                   opt.bc,...   
                                                opt.src, 'region3D', opt.region3D);    
   totmob = sum(mob, 2);
   %% made to add capillary pressure   
   if isfield(fluid,'pc'),
      pc=fluid.pc(state);
      gpc=zeros(size(G.cells.num));
      
      if isfield(fluid,'gpc') && strcmp(opt.pc_form,'global'),
         gpc=fluid.gpc(state);
      end

      if ~all(pc==0),
         if isfield(fluid,'gpc') && strcmp(opt.pc_form,'global'),
            cc = capPressureRHS(G,mob,pc,gpc,opt.pc_form); % ERROR HERE
         else
            %mob ok since we set bnd = 0 afterwards
            cc = capPressureRHS(G,mob,pc,opt.pc_form);
            % remove contribution on faces shared by 2D and 3D cells since
            % we use the 3D discretization here
            cc(G.facesBnd.cellFace2D) = 0;    
            cc(G.facesBnd.cellFace3D) = 0;  
         end
         grav = grav + cc;
      end
   end
   
   % gravity is defined on cellFaces, so we must acumulate to faces
      
   sgn = 2*(G.faces.neighbors(cf, 1) == cellNo) - 1;
   % cellFaces - rhs has contributions from internal and pressure faces
   j   = i(cf) | dF(cf);
   % on faces
   fg  = accumarray(cf(j), grav(j).*sgn(j), [nf, 1]);
   
   bndCells = reshape(G.faces.neighbors(bndFaces2D3D,:)',[],1);
   bndCells2Drep = bndCells(~opt.region3D(bndCells));
   sgn2D = 2*(G.faces.neighbors(bndFaces2D3D, 1) ==  bndCells2Drep)-1;
   
   % add gravity correction to faces, possibly wrong sign?!
   fg(bndFaces2D3D) = fg(bndFaces2D3D) + gravCorrVE.*sgn2D;

   % on cells
   rhs = accumarray(cellNo, -ft(cf).*(sgn.*fg(cf)+ff), [n, 1]) + ...
         [gg; zeros(nw, 1)]                                    + ...
         accumarray(cellNo, -hh(cf), [n, 1]);
      
   d  = zeros(G.cells.num, 1);

   % Wells --------------------------
   C    = cell (nw, 1);
   D    = zeros(nw, 1);
   W    = opt.wells;

   for k = 1 : nw,
      wc       = W(k).cells;
      nwc      = numel(wc);
      w        = k + nc;

      wi       = W(k).WI .* totmob(wc);

      dp       = norm(gravity()) * W(k).dZ*sum(rho .* W(k).compi, 2);
      d   (wc) = d   (wc) + wi;

      if     strcmpi(W(k).type, 'bhp'),
         ww=max(wi);
         %ww=1.0;
         rhs (w)  = rhs (w)  + ww*W(k).val; 
         rhs (wc) = rhs (wc) + wi.*(W(k).val + dp);
         C{k}     = -sparse(1, nc);
         D(k)     = ww;

      elseif strcmpi(W(k).type, 'rate'),
         rhs (w)  = rhs (w)  + W(k).val;
         rhs (wc) = rhs (wc) + wi.*dp;

         C{k}     =-sparse(ones(nwc, 1), wc, wi, 1, nc);
         D(k)     = sum(wi);

         rhs (w)  = rhs (w) - wi.'*dp;

      else
         error('Unsupported well type.');
      end
   end

   C = vertcat(C{:});
   D = spdiags(D, 0, nw, nw);
   %-----------------------------------------


   % Add up internal face transmissibilities plus Dirichlet pressure
   % faces for each cell.
   d  = d + ...
        accumarray(cellNo(dF(cf)), T(dF(cf)), [nc, 1]) +...
        accumarray(reshape(G.faces.neighbors(i,:), [], 1), ...
                   repmat(ft(i), [2,1]),  [nc, 1]);

   % Assemble coefficient matrix for internal faces.  Boundary conditions
   % may introduce additional diagonal entries.  Also, wells introduce
   % additional equations and unknowns.
   I  = [G.faces.neighbors(i,1); G.faces.neighbors(i,2); (1:nc)'];
   J  = [G.faces.neighbors(i,2); G.faces.neighbors(i,1); (1:nc)'];
   V  = [-ft(i); -ft(i); d];
   A  = sparse(double(I), double(J), V, nc, nc);
   A = [A, C'; C D];

   tocif(opt.Verbose);


   %% ---------------------------------------------------------------------
   dispif(opt.Verbose, 'Solving linear system...\t\t\t');
   ticif (opt.Verbose);

   if ~any(dF) && (isempty(W) || ~any(strcmpi({ W.type }, 'bhp'))),
      if A(1) > 0,
         A(1) = 2*A(1);
      else
         error('A(1) not > 0')
      end
   end
   if opt.condition_number,
      %tic;
      disp('***************************************');
      disp(['Conditon number is :   ', num2str(condest(A))]);
      disp('***************************************');
      %toc;
   end
   p = opt.LinSolve(A, rhs);

   tocif(opt.Verbose);

   %% ---------------------------------------------------------------------
   dispif(opt.Verbose, 'Computing fluxes, face pressures etc...\t\t');
   ticif (opt.Verbose);

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
   state.pressure(1 : nc) = p(1 : nc);
   state.flux(:)          = flux;
   state.facePressure     = fpress;

   for k = 1 : nw,
      wc       = W(k).cells;
      dp       = norm(gravity()) * W(k).dZ*sum(rho .* W(k).compi, 2);
      state.wellSol(k).flux     = W(k).WI.*totmob(wc).*(p(nc+k) + dp - p(wc));
      state.wellSol(k).bhp = p(nc + k);
   end

   if opt.MatrixOutput,
      state.A   = A;
      %state.rhs = rhs;
   end
   if isfield(fluid,'pc')
     % state.p_ph=calcPhasePressure(pc,gpc,opt.pc_form,state.pressure(1:nc));
   end
   tocif(opt.Verbose);
   
   
   % DEBUG
   state.mob = mob;
   state.mobCF = mobCF;
end

%--------------------------------------------------------------------------

function [gravCorrVE, mob, mobCF, omegaCF, rho] = dynamic_quantities(state, fluid, G, cellNo, opt)
   gravCorrVE = zeros(numel(G.facesBnd.index),1);
   
   [mu, rho] = fluid.properties(state);
   s         = fluid.saturation(state);
   kr        = fluid.relperm(s, state);

   mob    = bsxfun(@rdivide, kr, mu);
   
   mobCF = mob(cellNo, :);    
   
   if any(opt.region3D)
       % Find the cells which are neighbors to a 3D region (i.e. both 2D
       % and 3D faces)
       
       bndCells2D = G.facesBnd.cells2D;
       
       faceIx = mcolon(G.cells.facePos(bndCells2D),G.cells.facePos(bndCells2D+1)-1).';   

       %find z-coordinate of centroid of top face, i.e. where pressure is
       %defined
       topFace_z = G.faces.centroids(G.cells.faces(faceIx((G.cells.faces(faceIx,2)==5)),1), end);

       numSubCells = (G.cells.columnPos(G.cells.mapTopSurface(bndCells2D)+1)-...
       G.cells.columnPos(G.cells.mapTopSurface(bndCells2D)));

       cellNo2D = rldecode(bndCells2D, numSubCells,1);
       
       % distance from top surface to cell centroid
       dz_top = G.facesBnd.centroid3D(:, end)-rldecode(topFace_z, numSubCells, 1); 

       [h, h_max] = fluid.sat2height(state);
       
       if norm(gravity)>0
           % GRAVITY CORRECTION for pressure
           % Cells are filled to centroid with CO2
           isFilledToCent = dz_top < G.cells.z(cellNo2D) + h(cellNo2D);
           % To h, the mass above is density 1 (CO2). After that it will be
           % whatever density 2 is (second phase). If a cell is filled to
           % the centroid with CO2, use only CO2 for calculating gravity.
           gravCorrVE = h(cellNo2D).*rho(1) + max(0, dz_top-h(cellNo2D)).*rho(2);
           gravCorrVE(isFilledToCent) = dz_top(isFilledToCent).*rho(1);  
           % Shouldn't this possibly be just the z component?
           gravCorrVE = gravCorrVE*norm(gravity);

           %change from cellFaces to faces
           gravCorr = zeros(numel(G.faces.neighbors(:,1)),1);
           gravCorr(G.cells.faces(G.facesBnd.cellFace2D, 1)) = gravCorrVE;
           gravCorrVE = gravCorr(G.facesBnd.index); 
       end

       % Go via topgrid to get the height
      [s_full h_top] = normalizeValuesVE(G.topSurfaceGrid, state, fluid, 'CoupledGrid', G);
      [n, t] = fillDegree(h_top, G.topSurfaceGrid);
      fill2d = t(G.facesBnd.cells2D);
      
      n2D = numel(G.facesBnd.cells2D);
      Gt = G.topSurfaceGrid;
      cP = Gt.cells.columnPos;
      partial = false(numel(G.facesBnd.cellFace2D),1);
      filled  = partial;
      
      % Find partial, filled and empty cells.
      for i = 1:n2D
          ind = G.facesBnd.cells2D(i);
          cells = Gt.columns.cells(cP(ind):cP(ind+1)-1);
          partcell = G.facesBnd.map3D == cells(n(ind)+1);
          partial(partcell) = true;
          filled(sum(numSubCells(1:i-1)) + 1: find(partcell)) = true;
      end
      % Empty is true where there is no CO2
      % Filled is true where s_CO2 == 1
      % Partial is true where s_CO2 is not 0 or 1
      empty = ~filled;
      filled = filled & ~partial;
      fill = rldecode(fill2d, numSubCells);
      % Linear combination of saturations as a function of fill degree
      ph_fill = [fill (1-fill)];
      % See fixCoupledBoundaryMob for more details - this is essentially a
      % modification of half face mobilities to ensure that the transition
      % from 2D to 3D is okay since the saturation has different
      % interpretations on each side of the boundary.
      mobCF(G.facesBnd.cellFace2D(filled ),:) = repmat([1 0]./mu, sum(filled) ,1);
      mobCF(G.facesBnd.cellFace2D(empty  ),:) = repmat([0 1]./mu, sum(empty ) ,1);
      mobCF(G.facesBnd.cellFace2D(partial),:) = ph_fill(partial,:)./repmat(mu, sum(partial), 1);

   end   
   totmobCF = sum(mobCF, 2);
   omegaCF  = sum(bsxfun(@times, mobCF, rho), 2) ./ totmobCF;  
end
