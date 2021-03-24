function state = incompTPFA_DFM(state, G, T, fluid, varargin)
%Solve incompressible flow problem (fluxes/pressures) using TPFA method.
%
% SYNOPSIS:
%   state = incompTPFA_DFM(state, G, T, fluid)
%   state = incompTPFA_DFM(state, G, T, fluid, 'pn1', pv1, ...)
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
%   The function is modified from the original incompMPFA (legacy version) to account for
%   non-neighbor cell-to-cell connections for hybrid cells.
%
% NOTE:
%   This file is modified from the original incompTPFA to account for the
%   precence of hybrid cells in the grid. Furthermore, the mobility can be
%   treated by upstream weighting instead of harmonic averaging (only for
%   horizontal flow without capillary pressure).
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
%   cellConnection -A vector containing the transmissibilities for the cell
%                   connections. The connections must be spesified in
%                   G.cells.neighbors.
%                   A vector state.fluxc2c is returned contining the cell2cell
%                   fluxes.
%
%   upwind         -If true, the upwind mobility is used, els the harmonic
%                   averange of the mobility is used
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
%               - fluxc2c -- Fluxes over cell-to-cell connections
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
%   `computeTrans`, `addBC`, `addSource`, `addWell`, `initSingleFluid`, `initResSol`,
%   `initWellSol`.

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


   opt = struct('bc', [], 'src', [], 'wells', [], ...
                'LinSolve',     @mldivide,        ...
                'MatrixOutput', false, ...
                'Verbose',      mrstVerbose,...
                'condition_number',false,...
                'c2cTrans',[], ...
                'upwind',false,...
                'onlyMatrixOutput', false, ...
                'pc_form','wetting');

   opt = merge_options(opt, varargin{:});

   g_vec   = gravity();
   no_grav = ~(norm(g_vec(1 : size(G.nodes.coords,2))) > 0);
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

   % Preliminaries
   cellNo = rldecode(1:G.cells.num, double(diff(G.cells.facePos)), 2).';
   cf     = G.cells.faces(:,1);
   nf     = G.faces.num;
   nc     = G.cells.num;
   nw     = length(opt.wells);
   n      = nc + nw;

   % calculate the total mobility in each cell.
   [mob, omega, rho] = dynamic_quantities(state, fluid);
   totmob = sum(mob, 2);

   % Find cell-to-cell connections (presumably due to the elimination of
   % small cells in the intersection of fractures).
   if(~isempty(opt.c2cTrans))
       assert(isfield(G.cells,'neighbors'),'cell connections must be spesified')

       % The hybrid cell2cell transmissibilities
       ft_hyb = opt.c2cTrans;
   else
       ft_hyb = [];
   end

   % Identify internal faces
   i  = all(G.faces.neighbors ~= 0, 2);

   % Find the total transmissibility by multipling the mobility.
   % If specified the upwind mobility is used, else a harmonic average is used
   % for the mobilities.
   if opt.upwind
       % Face transmissibility = harmonic average of half-transmissibilities
       ft = 1 ./ accumarray(cf, 1./T, [nf, 1]);

       % pick out the upwind fluxes
       negFlow = state.flux < 0;
       N = G.faces.neighbors;
       upwind = zeros(length(N),1);
       upwind(~negFlow & i) = N(~negFlow & i,1);
       upwind(negFlow & i) = N(negFlow & i,2);
       upwind(~i) = sum(N(~i,:),2);

       % update the half and face transmissibility
       T  = T .* totmob(upwind(cf));
       ft = ft.*totmob(upwind);

       % find the upwind direction for the cell2cell connections.
       if ~isempty(opt.c2cTrans)
           negFlow = state.fluxc2c < 0;
           Nh = G.cells.neighbors;
           upwind2 = zeros(length(Nh),1);
           upwind2(~negFlow) = Nh(~negFlow,1);
           upwind2(negFlow)  = Nh(negFlow,2);

           % update the cell2cell face transmissibility
           ft_hyb = ft_hyb.*totmob(upwind2);
       end
   else
       % els the harmonic average of the total transmissibility is used.
       T = T .* totmob(cellNo);
       % Face transmissibility = harmonic average of half-transmissibilities
       ft = 1 ./ accumarray(cf, 1./T, [nf, 1]);

       % Star-delta transformation for the cell2cell connections
       if(~isempty(opt.c2cTrans))
           % this is an error start delta transformen
           % warning('Errors on large than 2 connections')
           [~,totmob_hyb] = computeHybridTrans(G,totmob(cellNo));
           ft_hyb 	= ft_hyb.*totmob_hyb;
       end
   end


   % Boundary conditions and source terms.
   [ff, gg, hh, grav, dF, dC] = computePressureRHS(G,...
                                                   omega, ...
                                                   opt.bc,...
                                                   opt.src);
   %% made to add capillary pressure
   if isfield(fluid,'pc'),
      pc=fluid.pc(state);
      gpc=zeros(size(totmob));

      if isfield(fluid,'gpc') && strcmp(opt.pc_form,'global'),
         gpc=fluid.gpc(state);
      end

      if ~all(pc==0),
         if isfield(fluid,'gpc') && strcmp(opt.pc_form,'global'),
            cc = capPressureRHS(G,mob,pc,gpc,opt.pc_form); % ERROR HERE
         else
            cc = capPressureRHS(G,mob,pc,opt.pc_form);
         end
         grav = grav + cc;
      end
   end

   sgn = 2*(G.faces.neighbors(cf, 1) == cellNo) - 1;
   j   = i(cf) | dF(cf);
   fg  = accumarray(cf(j), grav(j).*sgn(j), [nf, 1]);

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

   % Assamble the discretization matrix

   % If cell2cell connections are present additional equations and unknowns
   % are added.
   if(~isempty(opt.c2cTrans) && ~isempty(G.cells.neighbors))
       n1 = G.cells.neighbors(:,1);
       n2 = G.cells.neighbors(:,2);
       d = d + accumarray([n1;n2],[ft_hyb;ft_hyb],[nc,1]);
       I  = [G.faces.neighbors(i,1); G.faces.neighbors(i,2); n1; n2; (1:nc)'];
       J  = [G.faces.neighbors(i,2); G.faces.neighbors(i,1); n2; n1; (1:nc)'];
       V  = [-ft(i); -ft(i); -ft_hyb;-ft_hyb;d];
   else
       I  = [G.faces.neighbors(i,1); G.faces.neighbors(i,2); (1:nc)'];
       J  = [G.faces.neighbors(i,2); G.faces.neighbors(i,1); (1:nc)'];
       V  = [-ft(i); -ft(i);d];
   end

   A  = sparse(double(I), double(J), V, nc, nc);
   A = [A, C'; C D];



   tocif(opt.Verbose);
    if opt.MatrixOutput,
      state.A   = A;
      state.rhs = rhs;
   end

   if opt.onlyMatrixOutput, return, end

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

   % Dirichlet faces
   fpress(dF) = dC;

   % Sign for boundary faces
   sgn  = 2*(G.faces.neighbors(~i,2)==0)-1;
   ni   = G.faces.neighbors(i,:);
   flux = -accumarray(find(i),  ft(i) .*(p(ni(:,2))-p(ni(:,1))-fg(i)), [nf, 1]);
   c    = sum(G.faces.neighbors(~i,:),2) ;
   fg  = accumarray(cf, grav, [nf, 1]);
   flux(~i) = -sgn.*ft(~i).*( fpress(~i) - p(c) - fg(~i) );

   % Compute cell2cell flux
   if(~isempty(opt.c2cTrans))
       flux_hyb =-ft_hyb .*(p(n2)-p(n1));
   else
       flux_hyb=[];
   end

   % store the pressure and flux
   state.pressure(1 : nc) = p(1 : nc);
   state.flux(:)          = flux;
   state.fluxc2c = flux_hyb;
   state.facePressure     = fpress;


   for k = 1 : nw,
      wc       = W(k).cells;
      dp       = norm(gravity()) * W(k).dZ*sum(rho .* W(k).compi, 2);
      state.wellSol(k).flux     = W(k).WI.*totmob(wc).*(p(nc+k) + dp - p(wc));
      state.wellSol(k).pressure = p(nc + k);
   end

   if opt.MatrixOutput,
      state.A   = A;
      state.rhs = rhs;
   end
   if isfield(fluid,'pc')
      state.p_ph=calcPhasePressure(pc,gpc,opt.pc_form,state.pressure(1:nc));
   end
   tocif(opt.Verbose);
end

%--------------------------------------------------------------------------

function [mob, omega, rho] = dynamic_quantities(state, fluid)
   [mu, rho] = fluid.properties(state);
   s         = fluid.saturation(state);
   kr        = fluid.relperm(s, state);

   mob    = bsxfun(@rdivide, kr, mu);
   totmob = sum(mob, 2);
   omega  = sum(bsxfun(@times, mob, rho), 2) ./ totmob;
end
