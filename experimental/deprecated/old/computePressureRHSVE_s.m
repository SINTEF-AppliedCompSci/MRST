function [ff, gg, hh, gp, dF, dC] = computePressureRHSVE_s(g, omega, bc, src, varargin)
%Compute right-hand side contributions to pressure linear system for VE-s.
%
% SYNOPSIS:
%   [f, g, h, grav, dF, dC] = computePressureRHS_s(G, omega, bc, src)
%
% DESCRIPTION:
%   The contributions to the right-hand side for mimetic, two-point and
%   multi-point discretisations of the equations for pressure and total
%   velocity 
%                     
%     v + lam K·grad (p - g·z omega) = 0
%     div v  = q 
%
%   with  
%             __
%             \    kr_i/
%     lam   = /       / mu_i 
%             -- i
%             __
%             \
%     omega = /    f_i·rho_i
%             -- i
%              
% PARAMETERS:
%   G     - Grid data structure.
%
%   omega - Accumulated phase densities \rho_i weighted by fractional flow
%           functions f_i -- i.e., omega = \sum_i \rho_i f_i.  One scalar
%           value for each cell in the discretised reservoir model, G.
%
%   bc    - Boundary condition structure as defined by function 'addBC'.
%           This structure accounts for all external boundary conditions to
%           the reservoir flow.  May be empty (i.e., bc = struct([])) which
%           is interpreted as all external no-flow (homogeneous Neumann)
%           conditions.
%
%   src   - Explicit source contributions as defined by function
%           'addSource'.  May be empty (i.e., src = struct([])) which is
%           interpreted as a reservoir model without explicit sources.
%
% RETURNS:
%   f, g, h - Pressure (f), source/sink (g), and flux (h) external
%             conditions.  In a problem without effects of gravity, these
%             values may be passed directly on to linear system solvers
%             such as 'schurComplementSymm'.
%
%   grav    - Pressure contributions from gravity,
%
%                grav = omega·g·(x_face - x_cell)
%
%             where 
%
%                omega = \sum_i f_i\rho_i,
%
%             thus grav is a vector with one scalar value for each
%             half-face in the model (size(G.cells.faces,1)).
%
%   dF, dC  - Packed Dirichlet/pressure condition structure.  The faces in
%             'dF' and the values in 'dC'.  May be used to eliminate known
%             face pressures from the linear system before calling a system
%             solver (e.g., 'schurComplementSymm').
%
% SEE ALSO:
%   addBC, addSource, computeMimeticIP, schurComplementSymm.

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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

% $Date: 2012-01-30 11:41:03 +0100 (Mon, 30 Jan 2012) $
% $Revision: 9020 $    
   gp = grav_pressure(g, omega);
   ff = zeros(size(gp));
   gg = zeros([g.cells.num, 1]);
   hh = zeros([g.faces.num, 1]);

   if ~isempty(src),
      assert (max(src.cell) <= g.cells.num);

      ss = accumarray(src.cell, src.rate)    ;
      ii = accumarray(src.cell,    1    ) > 0;
      gg(ii) = gg(ii) + ss(ii);
   end

   dF = false([g.faces.num, 1]);
   dC = [];

   if ~isempty(bc),
      % Check that bc and g are compatible.
      assert (max(bc.face) <= g.faces.num && min(bc.face) > 0, ...
         'Boundary condition refer to face not existant in grid.');
      assert (all(accumarray(bc.face, 1, [g.faces.num, 1]) <= 1), ...
         'There are repeated faces in boundary condition.');

      % Pressure (Dirichlet) boundary conditions.
      %  1) Extract the faces marked as defining pressure conditions.
      %     Define a local numbering (map) of the face indices to the
      %     pressure condition values.
      %
      is_press = strcmpi('pressure', bc.type);
      face     = bc.face (is_press);
      dC       = bc.value(is_press);
      map      = sparse(double(face), 1, 1 : numel(face));

      %  2) For purpose of (mimetic) pressure solvers, mark the 'face's as
      %     having pressure boundary conditions.  This information will be
      %     used to eliminate known pressures from the resulting system of
      %     linear equations.  See (e.g.) 'solveIncompFlow' for details.
      %
      dF(face) = true;

      %  3) Enter Dirichlet conditions into system right hand side.
      %     Relies implictly on boundary faces being mentioned exactly once
      %     in g.cells.faces(:,1).
      %
      i     =   dF(    g.cells.faces(:,1) );
      ff(i) = - dC(map(g.cells.faces(i,1)));

      %  4) Reorder Dirichlet conditions according to SORT(face).  This
      %     allows the caller to issue statements such as 'X(dF) = dC' even
      %     when ISLOGICAL(dF).
      %
      dC = dC(map(dF));

      % Flux (Neumann) boundary conditions.
      % Note negative sign due to bc.value representing INJECTION flux.
      %
      is_flux              = strcmpi('flux', bc.type);
      hh(bc.face(is_flux)) = - bc.value(is_flux);
   end

   assert (~any(dC < 0));  % Pressure conditions should always be non-neg.
end

%--------------------------------------------------------------------------

function ff = grav_pressure(g, omega)
% computes innerproduct cf (face_centroid - cell_centroid) * g for each face

g_vec = gravity();   % Must be a 1-by-3 row vector for subsequent code.

   if norm(g_vec) > 0,       
      
      dim = size(g.nodes.coords,2);

      assert (1 < dim && dim < 4);
      assert (all(size(g_vec) == [1,3]));
      cellno = rldecode(1 : g.cells.num, diff(g.cells.facePos), 2) .';     
      
      if (any(strcmp(g.type, 'topSurfaceGrid')))
         % VE - 2D (i.e. g is created using topSurfaceGrid)
         cvec   = (g.faces.z(g.cells.faces(:,1), :) - ...
            g.cells.z(cellno            , :));
         ff     = omega(cellno) .* cvec*norm(gravity);
      else
         % VE - 3D
         cvec   = g.faces.centroids(g.cells.faces(:,1), :) - ...
                  g.cells.centroids(cellno            , :);
         cellF2D = (1:size(g.cells.faces,1))';
         
         
         hftb = any(bsxfun(@eq,g.cells.faces(cellF2D,2),[5,6]),2);
         % remove contribution from top and bottom, but only for 2D (VE)
         % cells i.e. not for real 3D cells (in the coupled version)         
         cvec(cellF2D(hftb),:) = 0;
         ff     = omega(cellno) .* (cvec * g_vec(1:dim).');
         
      end
   else
      ff     = zeros([size(g.cells.faces,1), 1]);
   end
end
