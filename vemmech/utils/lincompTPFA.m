function state = lincompTPFA(dt, state, G, T, pv, fluid, rock, varargin)
% Solve weakly compressible flow problem (fluxes/pressures) using TPFA method.
%
% SYNOPSIS:
%   state = lincompTPFA(state, G, T, fluid)
%   state = lincompTPFA(state, G, T, fluid, 'pn1', pv1, ...)
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
%            'lincompTPFA' and, possibly, a transport solver such as
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
%                  the caller of function 'lincompTPFA'.
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
%           'lincompTPFA:DrivingForce:Missing'
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
%    state = lincompTPFA(state, G, T, f, 'bc', bc, 'src', src, ...
%                       'wells', W, 'MatrixOutput', true);
%
%    plotCellData(G, state.pressure)
%
% SEE ALSO:
%   `computeTrans`, `addBC`, `addSource`, `addWell`, `initSingleFluid`, `initResSol`,
%   `initWellSol`.

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


   opt = struct('bc', [], 'src', [], 'wells', [], ...
                'bcp',[],...
                'LinSolve',     @mldivide,        ...
                'MatrixOutput', false, ...
                'Verbose',      mrstVerbose,...
                'condition_number',false,...
                'pc_form','nonwetting',...
                'use_trans',false);

   opt = merge_options(opt, varargin{:});

   g_vec   = gravity();
   % If gravity is overriden, we cannot say anything about the effects of
   % gravity on rhs.
   grav = (norm(g_vec(1:G.griddim)) > 0) || isfield(G, 'grav_pressure');

   if all([isempty(opt.bc)   , ...
           isempty(opt.src)  , ...
           isempty(opt.bcp)  , ...
           isempty(opt.wells), ~grav]),
      warning(msgid('DrivingForce:Missing'),                   ...
             ['No external driving forces present in model--', ...
              'state remains unchanged.\n']);
   end

   %% ---------------------------------------------------------------------
   dispif(opt.Verbose, 'Setting up linear system...\t\t\t');
   t0 = ticif (opt.Verbose);

   % Preliminaries
   [neighborship, n_isnnc] = getNeighbourship(G, 'Topological', true);
   [cellNo, cf, cn_isnnc] = getCellNoFaces(G);
   nif    = size(neighborship, 1);
   ncf    = size(cf, 1);
   nc     = G.cells.num;
   nw     = length(opt.wells);
   n      = nc + nw;

   [mob, omega, rho] = dynamic_quantities(state, fluid);
   totmob = sum(mob, 2);

   % Compute effective (mobility-weighted) transmissibilities.
   [T, ft] = compute_trans(G, T, cellNo, cf, neighborship, totmob, opt);

   % Identify internal faces
   i  = all(neighborship ~= 0, 2);

   % Boundary conditions and source terms.
   hh = zeros(nif, 1);
   dF = false(nif, 1);
   [grav, ff] = deal(zeros(ncf, 1));
   [ff(~cn_isnnc), gg, hh(~n_isnnc), grav(~cn_isnnc), dF(~n_isnnc), dC] = ...
      computePressureRHS(G, omega, opt.bc, opt.src);

   % made to add capillary pressure
   if isfield(fluid, 'pc'),
      pc  = fluid.pc(state);
      gpc = zeros(size(totmob));

      if isfield(fluid,'gpc') && strcmp(opt.pc_form,'global'),
         gpc = fluid.gpc(state);
      end

      if ~all(pc == 0),
         if isfield(fluid, 'gpc') && strcmp(opt.pc_form, 'global'),
            cc = capPressureRHS(G, mob, pc, gpc, opt.pc_form); % ERROR HERE
         else
            cc = capPressureRHS(G, mob, pc, opt.pc_form);
         end
         grav = grav + cc;
      end
   end
   clear mob

   sgn = 2*(neighborship(cf, 1) == cellNo) - 1;
   j   = i(cf) | dF(cf);
   fg  = accumarray(cf(j), grav(j).*sgn(j), [nif, 1]);
   if ~isempty(opt.bcp),
      fg(opt.bcp.face) = fg(opt.bcp.face) + opt.bcp.value;
      warning('mrst:periodic_bc', ...
             ['Face pressures are not well defined for ', ...
              'periodic boundary faces']);

      if any(G.faces.neighbors(:,1) == G.faces.neighbors(:,2)),
         error(['Periodic boundary: This code do not work ', ...
                'if a face is in and outflow']);
      end
   end


   rhs = accumarray(cellNo, -ft(cf).*(sgn.*fg(cf)+ff), [n, 1]) + ...
         [gg; zeros(nw, 1)]                                    + ...
         accumarray(cellNo, -hh(cf), [n, 1]); clear sgn;


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
        accumarray(reshape(neighborship(i,:), [], 1), ...
                   repmat(ft(i), [2,1]),  [nc, 1]);

   % Assemble coefficient matrix for internal faces.  Boundary conditions
   % may introduce additional diagonal entries.  Also, wells introduce
   % additional equations and unknowns.
   I  = [neighborship(i,1); neighborship(i,2); (1:nc)'];
   J  = [neighborship(i,2); neighborship(i,1); (1:nc)'];
   V  = [-ft(i); -ft(i); d]; clear d;
   A  = sparse(double(I), double(J), V, nc, nc);
   A = [A, C'; C D]; clear I J V C D;
   tocif(opt.Verbose, t0);


   %% ---------------------------------------------------------------------
   dispif(opt.Verbose, 'Solving linear system...\t\t\t');
   t0 = ticif (opt.Verbose);
%{
   if ~any(dF) && (isempty(W) || ~any(strcmpi({ W.type }, 'bhp'))),
      if A(1) > 0,
         A(1) = 2*A(1);
      else
         [j, j] = max(diag(A));  %#ok
         A(j,j) = 2 * A(j,j);
      end
   end
%}
   if opt.condition_number,
      disp('***************************************');
      disp(['Conditon number is :   ', num2str(condest(A))]);
      disp('***************************************');
   end
   %% add compressible effects
   nv=size(A,2);  
   state0=state;
   %pv=poreVolume(G,rock);
   if ~isfield(rock, 'cr')
       rock.cr = 0;
   end
   
   if ~isfield(fluid, 'cf')
       fluid.cf = 0;
   end
   cr=(rock.cr+fluid.cf);
   ct=cr.*pv/dt;
   rhs_c=ct.*state0.pressure;
   ct=sparse(1:G.cells.num,1:G.cells.num,ct,nv,nv);
   if opt.MatrixOutput,
      state.A   = A;
      state.rhs = rhs;
      state.ct = sparse(1:G.cells.num, 1 : G.cells.num, cr*pv, nv, nv);
   end
   A=A+ct;
   %orhs=rhs;
   rhs(1:nc)=rhs(1:nc)+rhs_c;
   %%
   
   p = opt.LinSolve(A, rhs);
   tocif(opt.Verbose, t0);

   clear A rhs;

   %% ---------------------------------------------------------------------
   dispif(opt.Verbose, 'Computing fluxes, face pressures etc...\t\t');
   t0 = ticif (opt.Verbose);

   % Reconstruct face pressures and fluxes.
   fpress     =  ...
          accumarray(cf, (p(cellNo)+grav).*T, [nif, 1])./ ...
          accumarray(cf(:,1), T, [nif,1]);


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
   ni   = neighborship(i,:);
   flux = -accumarray(find(i),  ft(i) .*(p(ni(:,2))-p(ni(:,1))-fg(i)), [nif, 1]);
   c    = sum(G.faces.neighbors(~i,:),2) ;
   fg  = accumarray(cf, grav, [nif, 1]);
   flux(~i) = -sgn.*ft(~i).*( fpress(~i) - p(c) - fg(~i) );
   %flux = -sgn.*ft((fpress(~i)-p(c)-grav));
   state.pressure(1 : nc) = p(1 : nc);
   state.flux(:)          = flux;
   state.facePressure     = fpress;

   for k = 1 : nw,
      wc       = W(k).cells;
      dp       = norm(gravity()) * W(k).dZ*sum(rho .* W(k).compi, 2);
      state.wellSol(k).flux     = W(k).WI.*totmob(wc).*(p(nc+k) + dp - p(wc));
      state.wellSol(k).pressure = p(nc + k);
   end

   tocif(opt.Verbose, t0);
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

%--------------------------------------------------------------------------

function [T, ft] = compute_trans(G, T, cellNo, cellFaces, neighborship, totmob, opt)
    niface = size(neighborship, 1);
    if opt.use_trans,  
      neighborcount = sum(neighborship > 0, 2);
      assert (numel(T) == niface, ...
             ['Expected one transmissibility for each interface ', ...
              '(=%d) but got %d.'], niface, numel(T));

      fmob = accumarray(cellFaces, totmob(cellNo), ...
                        [niface, 1]);
  
      fmob = fmob ./ neighborcount;
      ft   = T .* fmob;

      % Synthetic one-sided transmissibilities.
      th = ft .* neighborcount;
      T  = th(cellFaces(:,1));

   else

      % Define face transmissibility as harmonic average of mobility
      % weighted one-sided transmissibilities.
      %
      assert (numel(T) == numel(cellNo), ...
             ['Expected one one-sided transmissibility for each ', ...
              'half face (=%d), but got %d.'], numel(cellNo), numel(T));

      T  = T .* totmob(cellNo);
      ft = 1 ./ accumarray(cellFaces, 1 ./ T, [niface, 1]);

   end
end
