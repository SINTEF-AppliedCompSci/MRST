% Example setting up and solving a simple 2D poro-mechanics problem.
% See also mpfa_ex and mpsa_ex for more detailed comments.
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

% Construct grid
G = computeGeometry(cartGrid([2, 3], [1, 1]));

% Constitutive relation
mu = 1 : G.cells.num;
lambda = G.cells.num:-1:1;
constit = shear_normal_stress(G.cells.num, G.griddim, mu, lambda, 0 * mu);

% Boundary conditions
xmin_faces = [1, 4, 7];
xmax_faces = [3, 6, 9];
bc = addBC([], xmin_faces, 'pressure', 0); % Dirichlet on x=0

% Discretize elasticity problem and the Biot coupling terms
mpsa_discr = mpsa(G, constit, [], 'invertBlocks', 'matlab', 'bc', bc);

% boundary values
bc_vals = zeros(G.griddim * G.faces.num, 1); % Homogeneous BCs as default
% Face i has stress index 2*i-1, 2*i (2D)
bc_vals(2*xmin_faces) = 1:3; % Fix displacement in y-direction
bc_vals(2*xmax_faces-1) = -0.2; % Stress condition in x-direction

% Displacement is driven by boundary conditions
rhs_mech = -mpsa_discr.div * mpsa_discr.boundStress * bc_vals;


% Flow discretization
rock.perm = ones(G.cells.num, 1) * [10, 3, .5];% Perm: (kx=1, ky=3, kxy=.5)
mpfa_discr = mpfa(G, rock, [], 'invertBlocks', 'matlab', 'bc', bc);

% Some more parameters
compressibility = 1e-10; dt = 1; alpha_biot = 1;

% The Biot system is set up with elasticity in first row (block), flow in
% second
biot = [mpsa_discr.A , -alpha_biot * mpsa_discr.gradP ; ...
    alpha_biot * mpsa_discr.divD , ...
    alpha_biot^2 * mpsa_discr.stabDelta ... 
    + compressibility * spdiags(G.cells.volumes,0,G.cells.num,G.cells.num) ...
    + dt * mpfa_discr.A];

rhs = [rhs_mech; zeros(G.cells.num,1)]; % zero rhs for flow

u = biot \ rhs;

% The first Nd * num_cells numbers are displacements, these are again
% stored cell-wise
displacement = u(1:G.griddim * G.cells.num);
displacement_x = displacement(1:2:end);
displacement_y = displacement(2:2:end);
% Finally, the pressure
pressure = u(G.griddim * G.cells.num + 1 : end);