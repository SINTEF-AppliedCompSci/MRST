% Example setting up and solving a simple 2D elasticity problem using MPSA. 
% Much of the syntax is similar to that used for the flow problem
% (mpfa_ex.m), see there for more information.
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
% NOTE: If this does not work, it is most probably because MRST is not in
% your path.
G = computeGeometry(cartGrid([2, 3], [1, 1]));

% Constitutive relation
mu = 1 : G.cells.num;
lambda = G.cells.num:-1:1;
% shear_normal_stress assigns Lame-parameters to the data-structure for the
% constitutive law (one matlab-type cell per computational cell).
constit = shear_normal_stress(G.cells.num, G.griddim, mu, lambda, 0 * mu);

% Boundary conditions
xmin_faces = [1, 4, 7];
xmax_faces = [3, 6, 9];
bc = addBC([], xmin_faces, 'pressure', 0); % Dirichlet on x=0

% Compute discretization of elasticity. Biot coupling terms are also
% computed automatically, and will thus require some additional storage.
mpsa_discr = mpsa(G, constit, [], 'invertBlocks', 'matlab', 'bc', bc);

% Boundary values are assigned to a vector of size Nd * Num_cells. Face
% number one will have its x-condition assigned to bc_vals(1), y-condition
% at bc_vals(2). Next face no 2 etc.
bc_vals = zeros(G.griddim * G.faces.num, 1); % Homogeneous BCs as default
% Face i has stress index 2*i-1, 2*i (2D)
bc_vals(2*xmin_faces) = 1:3; % Fix displacement in y-direction
bc_vals(2*xmax_faces-1) = -0.2; % Stress condition in x-direction

% Displacement is driven by boundary conditions
rhs_mech = -mpsa_discr.div * mpsa_discr.boundStress * bc_vals;

% Solution and derived stress
displacement = mpsa_discr.A \ rhs_mech;
% Stress by cell center displacement and boundary conditions
stress = mpsa_discr.stress * displacement + mpsa_discr.boundStress * bc_vals;