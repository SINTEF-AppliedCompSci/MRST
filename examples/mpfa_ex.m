% Example setting up and solving a simple 2D flow problem using MPFA.
%
% The example assumes MRST is the Matlab path. For information on
% MRST-functions, confer the MRST documentation at
%   http://www.sintef.no/projectweb/mrst/
%
%{
Copyright 2015-2016, University of Bergen.

This file is part of FVBiot.

FVBiot is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

FVBiot is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%}


clear
% Initialize a Cartesian grid.
% NOTE: If this does not work, it is most probably because MRST is not in
% your path.
G = computeGeometry(cartGrid([20, 30], [1, 1]));
% Plot the grid using MRST functions
figure; plotGrid(G)

% Initialize permeability. Syntax and storage as a structure is chosen for
% compatibility with the equivalent method in MRST.
% rock.perm = ones(G.cells.num, 1) * [1, 3, .5];% Perm: (kx=1, ky=3, kxy=.5)
rock.perm = ones(G.cells.num, 1) * [1, 0, 1];% Perm: (kx=1, ky=3, kxy=.5)

% Indices of boundary faces

xmin_faces = find(G.faces.centroids(:, 1) == 0);
xmax_faces = find(G.faces.centroids(:, 1) == 1);


% xmin_faces = [1, 4, 7]; % Indices of faces on x=0
% xmax_faces = [3, 6, 9]; % Indices of faces at x=1

% Boundary conditions: Dirichlet on x=0, the rest will default to Neumann
% Note that the value given for the BC is somehow irrelevant (but needed
% for compatibility with MRST). We will instead assign the boundary values
% directly to bc_vals, below.
bc = addBC([], xmin_faces, 'pressure', 0);
bc = addBC(bc, xmax_faces, 'pressure', 0);

% The parameter invertBlocks describes which function is used for inverting
% local systems. This can be either 'matlab' or 'mex'. The mex option can
% give substantial speedup for large systems, in particular for elasticity,
% but it can also be memory-demanding.
mpfa_discr = mpfa(G, rock, [], 'invertBlocks', 'matlab', 'bc', bc);

% Assign boundary conditions
bc_vals = zeros(G.faces.num, 1);
bc_vals(xmin_faces) = 1/numel(xmin_faces)*(1 : numel(xmin_faces));
bc_vals(xmax_faces) = -0.2;

% Flow is driven by boundary conditions
rhs = -mpfa_discr.div * mpfa_discr.boundFlux * bc_vals;
pressure = mpfa_discr.A \ rhs;
% Plot pressure
plotCellData(G, pressure)
% Flux consists of contributions from pressure and bcs
flux = mpfa_discr.F * pressure + mpfa_discr.boundFlux * bc_vals;