function preComp = initTransportVE(G_top, rock2D, varargin)
% Precompute values needed in explicitTransportVE.
%
% SYNOPSIS:
%   preComp = initTransportVE(G_top, rock2D)
%
% PARAMETERS:
%   G_top  - structure representing the top-surface grid 3D grid.
%
%   rock2D - rock structure with porosities and (lateral) permeability
%            averaged for each column
%
% RETURNS:
%  preComp - structure of precomputed values containg the following fields:
%
%     - grav   - matrix of gravity flux contributions for each face/edge
%                for each cell. NB: weighted by 1/|c_ij| because we
%                multiply it by z_diff and h_diff to compute a term on the
%                form 'g*(grad z + grad h*rho_diff)'
%
%     - flux   - function handle for making matrix of flux contributions
%                for each face
%
%     - K_face - face permeability computed as a harmonic mean of the cell
%                permeabilites. Currently assumes K_x = K_y.
%
%     - pv     - pore volumes
%
%     - z_diff - vector of difference in z-coordinates for each for face ij
%                correspording to neighbors i,j:
%                z_diff(f_ij) = G.cells.z(i)-G.cells.z(j).
%
%     - n_z    - z component of unit normal of a cell. Used to
%                compute h_diff, perpendicular component.
%
%     - g_vec  - Vector of gravity flux for each face, used for computing
%                time steps and upwind mobility weighting
%
% COMMENTS:
%    Currently assumes rock.perm(:,1) = rock.perm(:,2)
%

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


% is_int - Logical vector of length G.faces.num, true for internal face.
is_int = all(double(G_top.faces.neighbors) > 0, 2);
[nc, nf] = deal(G_top.cells.num, sum(double(is_int)));
cellNo   = rldecode(1 : nc, diff(G_top.cells.facePos), 2) .';

% Indices to (internal) half-faces.
cIntFInx = is_int(G_top.cells.faces(:,1));

% Global-to-local face map for internal faces.
G2L      = cumsum(double(is_int));

% Indices of internal face corresponding to cIntFInx.
globfix  = G_top.cells.faces(cIntFInx, 1);

% Renumbering globfix to 1:numel(globfix)
locfix   = G2L(globfix);

magn = @(v)(sqrt(sum(v.^2,2)));

i = G_top.faces.neighbors(is_int,1);
j = G_top.faces.neighbors(is_int,2);

z_ij = G_top.cells.z(i)- G_top.cells.z(j);

% % vector between cell centroids
c_ij = [G_top.cells.centroids(i,:)-G_top.cells.centroids(j,:) z_ij];

% perpendicular component of unit normal of top surface
n_z    = G_top.cells.normals(:,3)./magn(G_top.cells.normals);

% face/edge-vector, all faces have 2 nodes: nb: must think about direction of e!
e_ij = zeros(G_top.faces.num, 3);
e_ij(:,1:2) = G_top.nodes.coords(G_top.faces.nodes(1:2:end-1),:) - ...
   G_top.nodes.coords(G_top.faces.nodes(2:2:end),:);
e_ij = e_ij(is_int,:);

g_cell = zeros([G_top.faces.num, 1]);
g_flux = abs(norm(gravity).*magn(e_ij)./magn(c_ij));
g_cell(is_int) = g_flux;

% - Gravity matrix -
% Spread gravity to cellfaces: cellFInx * g_const .* cellFace_normal
renum         = zeros([G_top.faces.num, 1]);
renum(is_int) = 1 : nf;

K2D = rock2D.perm(cellNo, 1);

% Compute harmonic average of nKg on all *internal* faces.
K_face = 2 ./ accumarray(renum(globfix), 1 ./ K2D(cIntFInx));

sgn = 2*double(G_top.faces.neighbors(G_top.cells.faces(:,1), 1) == cellNo) - 1;


% Compute final gravity matrix and make function for creating flux matrix
preComp.grav = sparse(cellNo(cIntFInx), locfix, ...
                      sgn(cIntFInx) .* g_cell(globfix), nc, nf);

preComp.flux = @(s) sparse(cellNo(cIntFInx), locfix, ...
                           sgn(cIntFInx) .*s.flux(globfix), nc, nf);
preComp.K_face = K_face;
preComp.pv     = rock2D.poro.*G_top.cells.volumes;
preComp.z_diff = z_ij;
preComp.n_z    = n_z;
preComp.g_vec  = g_flux;

end



