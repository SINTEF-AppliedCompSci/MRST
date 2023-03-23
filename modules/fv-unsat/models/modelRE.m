function [model] = modelRE(G, phys, mpfa_discr, bc, bcVal, relPermMethod, gEffects)
% Parent model for incompressible Richards' Equation 
%
% SYNOPSIS:
%   function [model] = modelRE(G, phys, mpfa_discr, bc, bcVal, relPermMethod, gEffects)
%
% PARAMETERS:
%   G             - Structure, Grid structure from MRST.
%   phys          - Structure, structure containing physical properties
%   krw_faces     - Vector, relative permeabilities at the faces
%   mpfa_discr    - Structure, mpfa discretization 
%   bc            - Structure, boundary conditions structure
%   bcVal         - Vector, containing the boundary condition values
%   relPermMethod - String, 'arithmetic' or 'upstream'
%   gEffects      - String, 'on' or 'off' to include/exclude gravity effects
%
% RETURNS:
%   model         - Structure, containing discrete equations and SWRC functions
%
% See also modelUnsatBiot.

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

% Importing modules
mrstModule add fv-unsat

% Grid-related quantities
zc = G.cells.centroids(:, end); % cell centers in z-direction
zf = G.faces.centroids(:, end); % face centers in z-direction
zetac = max(zf) - zc; % centroids of cells of elev. head
V = G.cells.volumes; % Cell volumes

% Gravity effects
if strcmp(gEffects, 'on')
    gravOn = 1;
elseif strcmp(gEffects, 'off')
    gravOn = 0;
else
    error('Gravity argument not recognized. Use either ''on'' or ''off''')
end

% Soil Water Retention Curves (SWRC) for the theta-psi model
[theta, krw, C_theta] = vGM_theta(phys);

% Discrete mpfa operators
F = @(x) mpfa_discr.F * x; % Flux  
boundF = @(x) mpfa_discr.boundFlux * x; % Boundary fluxes
divF = @(x) mpfa_discr.div * x; % Divergence

% Relative permeability at the faces
if strcmp(relPermMethod, 'arithmetic')
    krw_faces = @(psi_m) arithmeticAverageMPFA(G, krw, bc, psi_m);
elseif strcmp(relPermMethod, 'upstream')
    krw_faces = @(psi_m) upstreamWeightingMPFA(G, krw, bc, bcVal, ...
        mpfa_discr, phys, psi_m, 'psi', gEffects);
else
    error('Method not implemented. Use either ''arithmetic'' or ''upstream''')
end
        
% Darcy Flux
Q = @(psi, psi_m) (phys.flow.gamma ./ phys.flow.mu) .* krw_faces(psi_m) .* ...
    (F(psi + gravOn * zetac) + boundF(bcVal));

% Mass Conservation Equation                         
psiEq = @(psi, psi_n, psi_m, tau, source)  (V ./ tau) .* (theta(psi_m) ...
    + C_theta(psi_m) .* (psi - psi_m) - theta(psi_n)) ...
    + divF(Q(psi, psi_m)) - V .* source;

% Model equations
model = []; % initializing structure to store function handles
model.theta = theta; % Consitutive relationship for water content
model.krw = krw; % Consitutive relationship for rel perm
model.C_theta = C_theta; % Consitutive relationship specific moisture capacity
model.Q = Q; % Discrete Darcy equation
model.psiEq = psiEq; % Discrete Mass conservation