function pTop = computeTopPressure(G, phys, p, flux, modelEqs)
% Get mean pressure of the top layer using TPFA
%
% SYNOPSIS:
%   pTop = computeTopPressure(G, phys, p, flux, modelEqs)
%
% PARAMETERS:
%   G         - Structure, MRST grid structure
%   phys      - Structure, containing physical parameters
%   p         - Vector, containing the values of the pressure
%   flux      - Vector, Flux corresponding to the current m-level 
%   modelEqs  - Structure, containing the model equations
%
% RETURNS:
%   pTop      - Scalar, minumum value of top pressure
%

%{
Copyright 2018-2020, University of Bergen.

This file is part of the fv-unsat module.

fv-unsat is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

fv-unsat is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%} 

nz = G.numLayers; % number of layers
topCellsNum = G.cells.num / nz; % number of top cells
A  = G.faces.areas; % face areas
zc = G.cells.centroids(:, end); % cell centers in z-direction 
zf = G.faces.centroids(:, end); % face centers in z-direction
Lz = max(zf); % Depth of the domain
zetac = Lz - zc; % cell centers of elev. head
zetaf = Lz - zf; % face centers of elev. head
z_min = find(zf == 0); % idx of top faces

krw = modelEqs.krw; % relative permeability function
mu_w = phys.flow.mu; % viscosity
gamma = phys.flow.gamma; % specific gravity
k = mean(phys.flow.perm); % mean permeability

% Obtaining minimum value of pTop
pTop = mean( (((flux(z_min) ./ (A(z_min))) .* zc(1) .* mu_w) ./ ...
    (mean(krw(p(1:topCellsNum))) .* k)) - gamma .* ...
    (zetaf(z_min) - zetac(1:topCellsNum)) ...
    + p(1:topCellsNum));
end