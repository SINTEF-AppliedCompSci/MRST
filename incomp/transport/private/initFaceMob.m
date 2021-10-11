function [v_darcy, neighbors, ...
         normals, up_inx, g_vec] = initFaceMob(G, resSol, flux, grav)
%Initialize upwind face mobility saturation indices according to Darcy flux.
%That is, face mobility when not considering gravity.  Output is used by
%function findFaceMobMat.
%
% SYNOPSIS:
%   [v_darcy,   ...
%    neighbors, ...
%    normals,   ...
%    up_inx, g_vec] = initFaceMob(G, resSol, flux, grav)
%
%   G       - Grid data structure discretising the reservoir model.
%
%   resSol  - Reservoir solution structure containing valid (water)
%             saturation resSol.s(:,1) with one value for each cell in the
%             grid.
%
%   flux    - Inflow matrix of fluxes into each cell, created
%             by function initTransport.  Size
%             G.cells.num-by-(G.faces.num-#(boundary faces)). The fluxes
%             are assumed to be measured in units of m^3/s.
%
%   grav    - Matrix with gravity contribution for each face, created by
%             function initTransport.  Same size as flux.
%
% RETURNS:
%   v_darcy   - Darcy flux for internal faces.
%
%   neighbors - Neighbours for internal faces.
%
%   normals   - Vector of face normals for internal faces.
%
%   up_inx    - Index of cells for initial upwind face mobility.
%               Size = #(internal faces).
%
%   g_vec     - Gravity contributions = abs[K(rho_1-rho_2)*g*n_z],
%               for each face, where n_z is the z-component of the normal.
%
% SEE ALSO:
%   `findFaceMobMat`, `twophaseUpwFEGrav`, `twophaseUpwBEGrav`, `initTransport`.

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


intFInx   = all(double(G.faces.neighbors) > 0, 2);
neighbors = G.faces.neighbors(intFInx,:);
if nnz(flux),
   v_darcy = resSol.flux(intFInx);
else
   v_darcy = zeros([sum(double(intFInx)), 1]);
end

col    = 2 - double(~(v_darcy < 0));  % = 1 when v_darcy >= 0, 2 otherwise.
up_inx = sub2ind(size(G.faces.neighbors), find(intFInx), col);
up_inx = G.faces.neighbors(up_inx);

dim = size(G.faces.normals,2);
g   = gravity();

normals = sum(bsxfun(@times, g(1 : dim), G.faces.normals(intFInx,:)), 2);

use_grav = norm(g(1:dim)) > 0;
g_vec = 0;
if use_grav,
   % Compute vector of gravity contribution for each face.
   g_vec = sum(abs(grav), 1).' .* (0.5 * sign(normals));
end
