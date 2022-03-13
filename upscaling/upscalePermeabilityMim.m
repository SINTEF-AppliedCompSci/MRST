function perm_iris = ...
      upscalePermeabilityMim(Gp, bcp, dp_scale, Sp, fluid_pure, L)
% function for doing permeability upscaling on periodic grids
% perm_iris = upscalePermeability(Gp,G,bcp,dp_scale)
%
%  Gp       - period grid
%  bcp      - boundaryConditions
%  dp_scale - scale of pressure drop
%  Transp   - transmissibilities on periodic grid
%  fluid_pure - one phase fluid

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

   try
      require mimetic
   catch
      mrstModule add mimetic
   end

   dfaces = cell([Gp.griddim, 1]);
   vmat   = nan(Gp.griddim);

   for j = 1 : Gp.griddim,
      dfaces{j} = bcp.face(bcp.tags == j);
   end

   state   = initResSol(Gp, 100*barsa, 1);
   dp_mat  = dp_scale * eye(Gp.griddim);
   bcp_new = bcp;

   for i = 1 : Gp.griddim,
      for j = 1 : Gp.griddim,
         bcp_new.value(bcp.tags == j) = dp_mat(j,i);%/L(j);
      end

      state = incompMimetic(state, Gp, Sp, fluid_pure, 'bcp', bcp_new);

      for j = 1 : Gp.griddim,
         flux = sum(state.flux(dfaces{j}) .* bcp.sign(bcp.tags == j));

         vmat(j,i) = flux / sum(Gp.faces.areas(dfaces{j}));
      end
   end

   for j = 1 : Gp.griddim,
      dp_mat(j,:) = dp_mat(j,:) / L(j);
   end

   perm_iris = vmat / dp_mat;
end
