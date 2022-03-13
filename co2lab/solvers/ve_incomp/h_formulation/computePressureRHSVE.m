function [ff, gg, hh, gp, dF, dC] = computePressureRHSVE(g, omega, pc,  ...
                                                         bc, src, state)
%Compute right-hand side contributions to pressure linear system.
%
% SYNOPSIS:
%   [f, g, h, grav, dF, dC] = computePressureRHSVE(G, omega, bc, src, state)
%
% PARAMETERS:
%   G     - Grid data structure.
%
%   omega - Accumulated phase densities \rho_i weighted by fractional flow
%           functions f_i -- i.e., omega = \sum_i \rho_i f_i.  One scalar
%           value for each cell in the discretised reservoir model, G.
%
%   pc    - second-order term in the pressure equation ("capillary
%           pressure" function, fluid.pc).
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
%   state - Reservoir solution structure holding reservoir state.
%
% RETURNS:
%   f, g, h - Pressure (f), source/sink (g), and flux (h) external
%             conditions.  In a problem without effects of gravity, these
%             values may be passed directly on to linear system solvers
%             such as 'schurComplementSymm'.
%
%   grav    - Pressure contributions from gravity.  One scalar value for
%             each half-face in the model (size(G.cells.faces,1)).
%
%   dF, dC  - Packed Dirichlet/pressure condition structure.  The faces in
%             'dF' and the values in 'dC'.  May be used to eliminate known
%             face pressures from the linear system before calling a system
%             solver (e.g., 'schurComplementSymm').
%
% SEE ALSO:
%   `addBC`, `addSource`, `computeMimeticIP`, `schurComplementSymm`.

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

% $Date: 2012-01-30 11:41:03 +0100 (Mon, 30 Jan 2012) $
% $Revision: 9020 $

   gp = grav_pressure(g, omega, pc, state,bc);
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
      assert (max(bc.face) <= g.faces.num);
      assert (all(accumarray(bc.face, 1, [g.faces.num, 1]) <= 1));

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

%   assert (~any(dC < 0));  % Pressure conditions should always be non-neg.
end

%--------------------------------------------------------------------------

function ff = grav_pressure(g, omega, pc, state,bc)
% express (grad z g_parallel + grad_perp pc(h) g_perp) as a source term on
% each cellface

   if norm(gravity) > 0,
      dim   = size(g.nodes.coords,2);

      % unit cell normal
      cellnormals = bsxfun(@rdivide,g.cells.normals, ...
                           sqrt(sum(g.cells.normals.^2,2)));
      phi_cells = (g.cells.z+pc(state.h).*cellnormals(:,3))*norm(gravity());

      phi_hface = rldecode(phi_cells,diff(g.cells.facePos));
      % face value is average of value in cell1 cell2
      phi_face  = accumarray(g.cells.faces(:,1),phi_hface) ...
                  ./sum(g.faces.neighbors>0,2);

      %-- special treatment of external faces --%
      external      = sum(g.faces.neighbors>0,2)==1;
      external_cell = sum(g.faces.neighbors(external,:),2);

      % no flow bc: assume height is the cell height if closed boundary
      tmp_pc=pc(state.h);
      phi_face(external) = (g.faces.z(external)- tmp_pc(external_cell) ...
                           .*cellnormals(external_cell,3))*norm(gravity());
      % faces with other bc than no-flow
      if(~isempty(bc))
         bc_cell = sum(g.faces.neighbors(bc.face,:),2);
         phi_face(bc.face) =...
            (g.faces.z(bc.face)-pc(bc.h).*cellnormals(bc_cell,3))*norm(gravity());
      end

      assert(dim==2);

      cellno = rldecode(1 : g.cells.num, diff(g.cells.facePos), 2) .';

      ff     = omega(cellno) .* (phi_face(g.cells.faces(:,1))-phi_cells(cellno));
      %.*g.faces.areas(g.cells.faces(:,1));
   else
      ff     = zeros([size(g.cells.faces,1), 1]);
   end
end
