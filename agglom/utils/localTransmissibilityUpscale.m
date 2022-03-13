function T_coarse = localTransmissibilityUpscale(CG, trans, varargin)
%Use local transmissibility upscaling to find non-negative transmissibilities
%
% SYNOPSIS:
%   T_coarse = localTransmissibilityUpscale(CG, trans)
%   T_coarse = localTransmissibilityUpscale(CG, trans, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function finds a set of effective transmissibilities for each
%   coarse interface. By using a TPFA discretization and reasonable
%   boundary conditions, it *should* guarantee non-negative
%   transmissibilities for all coarse grids.
%
%   The algorithm used is based on attempting to make normal surfaces
%   parallel to the interface at each side of a pair of coarse blocks.
%   Cells further away than these planes are then assigned either pressure
%   1 or 0, which is then solved locally with no-flow elsewhere. The
%   resulting pressure is used to estimate the effective
%   transmissibilities.
%
% REQUIRED PARAMETERS:
%   CG    - Coarsegrid with valid geometry fields from coarsenGeometry.
%
%   trans - Prescribed fine-scale transmissibilities.  One positive scalar
%           for each connection.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   targetCoarseFaces -
%              Solve only for a subset of the coarse interfaces.  Note that
%              this number includes boundary faces on the coarse grid.
%              Either a logical mask of length CG.faces.num or a set of
%              indices in the range 1:CG.faces.num.  Default to all faces.
%
%   strict   - Strict will assert that all transmissibilities are
%              non-negative.  If this is disabled, the routine will not
%              check.  Default: ON.
%
%   returnInnerTrans -
%              Return only transmissibility for inner coarse faces (i.e.
%              all(CG.faces.neighbors > 0, 2)). Default: ON.
%
%   kinkless - The algorithm sets up two planes for each coarse face to
%              find boundary conditions.  The default is that kinkless is
%              off, resulting in each normal vector being defined from one
%              of the coarse block centroids to the coarse interface.  The
%              alternative is to set the normal vector from one block
%              centroid to the other, which is less robust, but may better
%              approximate a linear pressure drop. Default: FALSE.
%
% RETURNS:
%   T - Transmissibilities per inner coarse face (if returnInnerTrans is on)
%       or otherwise per global coarse face.  If set per global interface,
%       it will be zero at the boundary.

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

   internalFaces = all(CG.faces.neighbors > 0, 2);

   opt = struct('targetCoarseFaces', find(internalFaces), ...
                'strict'           , true               , ...
                'returnInnerTrans' , true               , ...
                'kinkless'         , false);

   opt = merge_options(opt, varargin{:});

   % Ensure that we have a row vector that we can iterate on
   opt.targetCoarseFaces = reshape(opt.targetCoarseFaces, 1, []);

   if islogical(opt.targetCoarseFaces)
      % Check for a logical mask of correct length
      assert(numel(opt.targetCoarseFaces) == CG.faces.num);
      opt.targetCoarseFaces = find(opt.targetCoarseFaces);
   end

   % Set up a incompressible equation system for the whole domain with noflow
   % at boundary.
   G = CG.parent;
   state = initResSol(G, 0, 1);

   % Unit fluid, does not really matter what the viscosity is.
   fluid = initSingleFluid('mu', 1, 'rho', 1);

   % Toggle off gravity and set up the system via incompTPFA
   oldGrav = gravity();
   gravity off

   warning('off', 'incompTPFA:DrivingForce:Missing')
   state = incompTPFA(state, G, trans, fluid, 'MatrixOutput', true, ...
                      'use_trans', true);
   warning('on', 'incompTPFA:DrivingForce:Missing')

   gravity(oldGrav)
   if norm(oldGrav) > 0,
      gravity on
   end

   A = state.A;

   % Set up face transmissibilities whilst incorporating non-neighbouring
   % connections (NNC).
   N = double(G.faces.neighbors);

   % Preallocate the internal faces
   T_coarse = zeros([CG.faces.num, 1]);

   % Make the coarsegrid lookups prettier
   cpos = CG.faces.connPos;
   fconn = CG.faces.fconn;
   for cf = opt.targetCoarseFaces,
      n = CG.faces.neighbors(cf, :);

      if any(n == 0),
         % Skip any boundary coarse faces
         continue
      end

      % Cells we are interested in belong to either "left" or "right"
      % side of the coarse interface
      pick = false([CG.cells.num, 1]); pick(n) = true;
      potentialCells = pick(CG.partition);  clear pick

      % Grab all centroids for categorization
      points = G.cells.centroids(potentialCells, :);

      % Take the centroids and the interface centroid to define planes
      fcent = CG.faces.centroids(cf);
      c1 = CG.cells.centroids(n(1), :);
      c2 = CG.cells.centroids(n(2), :);

      % Define normal vectors for planes, *always* oriented towards the
      % interface
      if opt.kinkless
         % We do not want a kink in our flow pattern, ignore the face
         % centroid.
         v1 = c2 - c1;
         v2 = -v1;
      else
         % Each side defined by itself
         v1 = fcent - c1;
         v2 = fcent - c2;
      end

      % We are interested in the cells between the two planes
      sgnA = checkPlaneSign(points, v1, c1) > 0;
      sgnB = checkPlaneSign(points, v2, c2) > 0;

      % Store both a logical and a actual lookup to avoid costly
      % operations later on
      fpcells = find(potentialCells);
      % Global lookup of internal nodes
      internal = fpcells(sgnA & sgnB);

      if ~ isempty(internal),
         % Exploit the incompressible global system to create a local
         % linear system for the internal nodes, with remaining
         % connections to all potential nodes and noflow to the rest of
         % the domain.
         A_local = A(internal, internal);
         A_local = A_local + diag(sum(A(internal, ~potentialCells), 2));

         % Define boundary condition by arbitrarily setting the outside of
         % the first plane to 1.  The others are implicitly set to zero
         % because we did not eliminate those equations.
         bnd = false(G.cells.num, 1);

         % Set boundary to those cells outside of plane 1
         bnd(fpcells(~sgnA)) = true;

         % But do not set boundary for cells belonging to the coarse block
         % not affiliated with plane 1.
         bnd(CG.partition == n(2)) = false;

         % rhs = 0 - A_(interface, bnd) * p_bnd =
         % A_(interface, leftbnd)*1...
         rhs = -A(internal, bnd) * ones([sum(bnd), 1]);

         % Solve problem and mask into global array because of limited time
         p_global           = zeros([G.cells.num, 1]);
         p_global(internal) = A_local \ rhs;

         % Set boundary condition in the strange case that those regions are
         % to be included in the interface (very small or thin coarse
         % blocks).
         p_global(bnd) = 1;
      else
         % If no cells are flagged as internal, it is likely that we
         % have blocks inside blocks or other highly degenerate cases.
         % In this case we fall back to using the sum of
         % transmissibilities over the interface. This is done by
         % setting p = 1 in one block and 0 in the other.
         p_global = double(ismember(CG.partition, n(1)));
         warning(['Possible block in block or other degenerate case ', ...
                  'for interface %d in coarse block pair (%d, %d).\n', ...
                  ' Falling back to sum of transmissibilities.', ...
                  ' Consider running ''removeConfinedBlocks'''], ...
                 cf, n(1), n(2));
      end

      % Work out fluxes over the interface
      conn = fconn(cpos(cf) : cpos(cf + 1) - 1);
      connTrans = trans(conn);

      % Find fine neighbors to each fine face in the interface
      n_iface = N(conn, :);

      % Set sign based on membership
      sgn = -1 + 2*ismember(n_iface(:, 1), find(CG.partition == n(1)));
      p1 = p_global(n_iface(:, 1));
      p2 = p_global(n_iface(:, 2));

      % The transmissibility is then the pressure difference times the
      % transmissibilities, summed over the whole coarse interface
      T_coarse(cf) = sum(connTrans .* (p1 - p2) .* sgn);
   end

   if opt.returnInnerTrans
      % Mask away boundary faces and return only for the inner interfaces
      T_coarse = T_coarse(internalFaces);
   end

   assert (~any(T_coarse < 0) || ~opt.strict, ...
           'Negative transmissibilities found, aborting!')
end

% --------------------------------------------------------------------------

function sgn = checkPlaneSign(points, normal, centroid)
% Given a set of points and a plane defined by a centroid and a normal
% vector, this defines the sign of any points with respect to that
% plane
   ptvec = bsxfun(@minus, points, centroid);
   sgn = sign(sum(bsxfun(@times, ptvec, normal), 2));
end
