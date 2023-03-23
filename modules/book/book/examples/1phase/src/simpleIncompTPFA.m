function state = simpleIncompTPFA(state, G, hT, fluid, varargin)
%Solve incompressible flow problem (fluxes/pressures) using TPFA method.
%
% SYNOPSIS:
%   state = simpleIncompTPFA(state, G, hT, fluid)
%   state = simpleIncompTPFA(state, G, hT, fluid, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function assembles and solves a (block) system of linear equations
%   defining interface fluxes and cell pressures at the next time step in a
%   sequential splitting scheme for the reservoir simulation problem
%   defined by Darcy's law and a given set of external influences (sources,
%   and boundary conditions).
%
%   This function uses a two-point flux approximation (TPFA) method with
%   minimal memory consumption within the constraints of operating on a
%   fully unstructured polyhedral grid structure.
%
% REQUIRED PARAMETERS:
%   state  - Reservoir solution structure either properly
%            initialized from function 'initResSol'
%
%   G, hT  - Grid and half-transmissibilities as computed by the function
%            'computeTrans'.
%
%   fluid  - Fluid object as defined by function 'initSimpleFluid'.
%
% OPTIONAL PARAMETERS:
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
%   LinSolve - Handle to linear system solver software to which the
%            fully assembled system of linear equations will be passed.
%            Assumed to support the syntax
%
%              x = LinSolve(A, b)
%
%            in order to solve a system Ax=b of linear equations.
%            Default value: LinSolve = @mldivide (backslash).
%
% RETURNS:
%   state - Update reservoir solution structure with new values
%           for the fields:
%              - pressure -- Pressure values for all cells in the
%                            discretised reservoir model, 'G'.
%              - facePressure --
%                            Pressure values for all interfaces in the
%                            discretised reservoir model, 'G'.
%              - flux     -- Flux across global interfaces corresponding to
%                            the rows of 'G.faces.neighbors'.
%
% NOTE:
%   If there are no external influences, i.e., if all of the structures
%   'bc' and 'src' are empty and there are no effects of gravity, then the
%   input value 'xr' is returned unchanged and a warning is printed in the
%   command window. This warning is printed with message ID
%
%           'incompTPFA:DrivingForce:Missing'
%
%
% SEE ALSO:
%   computeTrans, addBC, addSource, addWell, initSingleFluid, initResSol,
%   initWellSol.

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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

   opt = struct('bc', [], 'src', [], ...
                'LinSolve',     @mldivide);
   opt = merge_options(opt, varargin{:});

   % Sanity checks on grid and gravity term
   assert (1 <= G.griddim && G.griddim < 4);
   gvec = gravity(); gvec = gvec(1:G.griddim);
   assert (all(size(gvec) == [1,G.griddim]));
   if all([isempty(opt.bc)   , ...
           isempty(opt.src)  , ~(norm(gvec) > 0)]),
      warning(msgid('DrivingForce:Missing'),                   ...
             ['No external driving forces present in model--', ...
              'state remains unchanged.\n']);
   end

   % Preliminaries: set various constants and maps
   nf     = G.faces.num;
   nc     = G.cells.num;
   cf     = G.cells.faces(:,1);
   nhf    = numel(cf);
   hf2cn  = gridCellNo(G);
   iface  = all(G.faces.neighbors ~= 0, 2);
   eface  = ~iface;
   
   % Define effective face transmissibility as harmonic average of
   % viscosity-weighted one-sided transmissibilities.
   [mu, rho] = fluid.properties(state);
   assert (numel(hT) == numel(hf2cn), ...
      ['Expected one one-sided transmissibility for each ', ...
      'half face (=%d), but got %d.'], numel(hf2cn), numel(hT));
   hT  = hT ./ mu;
   T = 1 ./ accumarray(cf, 1 ./ hT, [G.faces.num, 1]);

   % Compute gravity contribution to right-hand side
   cvec = G.faces.centroids(cf, :) - G.cells.centroids(hf2cn, :);
   gp   = rho .* (cvec * gvec.'); clear cvec;
   
   % Initiatlize right-hand side
   rhs1 = zeros(nhf, 1);
   rhs2 = zeros(nc,  1);
   rhs3 = zeros(nf,  1);
   
   % Source terms
   src = opt.src;
   if ~isempty(src),
      % Compatibility check on cell numbers for source terms
      assert (max(src.cell) <= nc && min(src.cell>0), ...
         'Source terms refer to cell not existant in grid.');

      % Sum source terms inside each cell and add to rhs
      ss = accumarray(src.cell, src.rate);
      ii = accumarray(src.cell, 1) > 0;
      rhs2(ii) = rhs2(ii) + ss(ii);
   end

   % Boundary conditions
   Dface    = false([nf, 1]);
   DfaceVal = [];
   bc       = opt.bc;
   if ~isempty(bc),
       % Compatibility checks
      assert (max(bc.face) <= nf && min(bc.face) > 0, ...
         'Boundary condition refer to face not existant in grid.');
      assert (all(accumarray(bc.face, 1, [nf, 1]) <= 1), ...
         'There are repeated faces in boundary condition.');
      
      % Pressure (Dirichlet) boundary conditions.
      % Extract the faces marked as defining pressure conditions. Define a
      % local numbering (map) of the face indices to the pressure condition
      % values.
      is_press = strcmpi('pressure', bc.type);
      pface    = bc.face (is_press);
      DfaceVal = bc.value(is_press);
      map      = sparse(double(pface), 1, 1 : numel(pface));
      
      % Mark the faces as having pressure boundary conditions.  This
      % information will be used to eliminate known pressures from the
      % resulting system of linear equations.
      Dface(pface) = true;

      % Enter Dirichlet conditions into system right hand side. Relies
      % implictly on boundary faces being mentioned exactly once in
      % G.cells.faces(:,1).
      ind  = Dface(cf);
      rhs1(ind) = - DfaceVal(map(G.cells.faces(ind,1))); clear ind

      % Reorder Dirichlet conditions according to SORT(pface) so that we
      % later can set 'X(Dface) = DfaceVal' even when ISLOGICAL(Dface).
      DfaceVal = DfaceVal(map(Dface));

      % Flux (Neumann) boundary conditions.
      % Note negative sign due to bc.value representing INJECTION flux.
      is_flux = strcmpi('flux', bc.type);
      rhs3(bc.face(is_flux)) = - bc.value(is_flux);
   end
   assert(~any(DfaceVal < 0), 'Pressure conditions should always be non-negative');
   
   % Add gravity contribution to all internal faces and faces with
   % Dirichlet boundary conditions
   sgn = 2*(G.faces.neighbors(cf, 1) == hf2cn) - 1;
   j   = iface(cf) | Dface(cf);
   faceGrav  = accumarray(cf(j), gp(j).*sgn(j), [nf, 1]); clear j

   rhs = accumarray(hf2cn, -T(cf).*(sgn.*faceGrav(cf)+rhs1), [nc, 1]) + ...
      rhs2 + accumarray(hf2cn, -rhs3(cf), [nc, 1]); 
   clear rhs1 rhs2 sgn;


   % Add up internal face transmissibilities plus Dirichlet pressure
   % faces for each cell.
   d  = accumarray(hf2cn(Dface(cf)), hT(Dface(cf)), [nc, 1]) +...
        accumarray(reshape(G.faces.neighbors(iface,:), [], 1), ...
                   repmat(T(iface), [2,1]),  [nc, 1]);

   % Assemble coefficient matrix for internal faces.  Boundary conditions
   % may introduce additional diagonal entries. 
   I  = [G.faces.neighbors(iface,1); G.faces.neighbors(iface,2); (1:nc)'];
   J  = [G.faces.neighbors(iface,2); G.faces.neighbors(iface,1); (1:nc)'];
   V  = [-T(iface); -T(iface); d]; clear d;
   A  = sparse(double(I), double(J), V, nc, nc); clear I J V;

   % If there are no Dirichlet boundary conditions, do a fix to ensure that
   % we have a solvable system
   if ~any(Dface)
      if A(1) > 0,
         A(1) = 2*A(1);
      else
         [j, j] = max(diag(A));  %#ok
         A(j,j) = 2 * A(j,j);
      end
   end
   
   % Solve the flow problem
   p = opt.LinSolve(A, rhs);
   clear A rhs;

   % Reconstruct face pressures and fluxes.
   fpress =  ...
      accumarray(cf, (p(hf2cn)+gp).*hT, [G.faces.num,1])./ ...
      accumarray(cf, hT, [G.faces.num,1]);

   % Recompute face pressure at Neumann faces
   b         = any(G.faces.neighbors==0, 2);
   fpress(b) = fpress(b) - rhs3(b)./T(b);

   % Reset correct values at Dirichlet faces
   fpress(Dface) = DfaceVal;

   % Compute face fluxes for internal faces
   ni   = G.faces.neighbors(iface,:);
   flux = -accumarray(find(iface),  ...
      T(iface) .*(p(ni(:,2))-p(ni(:,1))-faceGrav(iface)), [nf, 1]);

   % Compute fluxes for external faces using Darcy's law
   sgn         = 2*(G.faces.neighbors(eface,2)==0)-1;
   cNo         = sum(G.faces.neighbors(eface,:),2) ;
   faceGrav    = accumarray(cf, gp, [nf, 1]);
   flux(eface) = -sgn.*T(eface).*( fpress(eface) - p(cNo) - faceGrav(eface) );
   
   % Extract the corresponding state variables
   state.pressure(1:nc) = p(1:nc);
   state.flux(:)        = flux;
   state.facePressure   = fpress;
end
