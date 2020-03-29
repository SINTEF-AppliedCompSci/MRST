function [model] = modelUnsatBiot(G, phys, mpfa_discr, mpsa_discr, ...
    bcFlow, bcFlowVals, bcMech, bcMechVals, relPermMethod, gEffects)
% Parent model for the equations of unsaturated poroelasticity
%
% SYNOPSIS
%   function [model] = modelUnsatBiot(G, phys, mpfa_discr, mpsa_discr, ...
%       bcFlow, bcFlowVals, bcMech, bcMechVals, relPermMethod, gEffects)
%
% PARAMETERS:
%   G             - Structure, Grid structure from MRST.
%   phys          - Structure, structure containing physical properties
%   mpfa_discr    - Structure, mpfa discretization 
%   mpsa_discr    - Structure, mpsa discretization 
%   bcFlow        - Structure, flow boundary condition
%   bcFlowVals    - Vector, containing flow boundary condition values
%   bcMech        - Structure, mechanics boundary condition
%   bcMechVals    - Vector, containing mechanics boundary condition values
%   relPermMethod - String, 'arithmetic' or 'upstream'
%   gEffects      - String, 'on' or 'off' to include/exclude gravity effects
%
% RETURNS:
%   model         - Structure, containing discrete equations and SWRC functions
%
% See also modelRE.

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


mrstModule add fv-unsat

% Grid-related quantities
Nc = G.cells.num; % total number of cells
Nd = G.griddim; % grid dimension
V  = G.cells.volumes; % cell volumes
zc = G.cells.centroids(:, end); % cell centers in z-direction     
zf = G.faces.centroids(:, end); % face centers in z-direction
zetac = max(zf) - zc; % cell centers of elev. head

% Gravity effects
if strcmp(gEffects, 'on')
    gravOn = 1;
elseif strcmp(gEffects, 'off')
    gravOn = 0;
else
    error('Gravity argument not recognized. Use either ''on'' or ''off''')
end

% Physical properties
flow = phys.flow;
mech = phys.mech;

% Soil Water Retention Curves (SWRC)
[S_w, krw, C_S] = vGM_saturation(phys);

% MPSA operators
S = @(x) mpsa_discr.stress * x; % stress
boundS = @(x) mpsa_discr.boundStress * x; % stress boundary
gradP = @(x) mpsa_discr.gradP * x; % gradient of pressure
divU = @(x) mpsa_discr.divD * x; % divergence of displacement
compat = @(x) mpsa_discr.stabDelta * x; % stability operator
divS = @(x) mpsa_discr.div * x; % divergence of stress

% MPFA operators
F = @(x) mpfa_discr.F * x; % flux
boundF = @(x) mpfa_discr.boundFlux * x; % boundary fluxes
divF = @(x) mpfa_discr.div * x; % divergence of flux

% Relative permeability at the faces
if strcmp(relPermMethod, 'arithmetic')
    krw_faces = @(p_m) arithmeticAverageMPFA(G, krw, bcFlow, p_m);
elseif strcmp(relPermMethod, 'upstream')
    krw_faces = @(p_m) upstreamWeightingMPFA(G, krw, bcFlow, bcFlowVals, ...
        mpfa_discr, phys, p_m, 'pressure', gEffects);
else
    error('Method not implemented. Use either ''arithmetic'' or ''upstream''')
end

% Mechanical discrete equations

% Function that maps "c" to "Nd*c"
sca2vec = @(x) repmat(x, [Nd, 1]); 

% Body forces
g_body = zeros(Nd * Nc, 1);  % intializing body forces vector
g_body(Nd:Nd:end) = flow.g;  % assigning g to z-cells
body = @(p_n) ((1 - flow.poro) .* mech.rho + ...
    flow.poro .* sca2vec(S_w(p_n)) .* flow.rho) .* g_body * gravOn;

% Traction
T = @(u) S(u) + boundS(bcMechVals);

% Momentum equation (Mechanics contribution)
uEq1 = @(u) divS(T(u));

% Momentum equation (Flow contribution + source)
uEq2 = @(p, p_n, sourceMech) -flow.alpha_biot .* gradP(S_w(p_n) .* p) ...
    + sca2vec(V) .* sourceMech;

% Flow discrete equations

% Compressibility-like terms
xi  = @(p_n) (flow.alpha_biot - flow.poro) .* mech.C_s .* S_w(p_n).^2 ...
    + flow.poro .* flow.C_w .* S_w(p_n);

chi = @(p_n) (flow.alpha_biot - flow.poro) .* mech.C_s .* S_w(p_n) .* p_n ...
    + flow.poro;  

% Darcy Flux
Q = @(p, p_m) (krw_faces(p_m) ./ flow.mu) .* ...
    (F(p + gravOn .* flow.gamma .* zetac) + boundF(bcFlowVals)); 

% Mass conservation (Mechanics contribution)
pEq1 = @(p_n, u, u_n) flow.alpha_biot .* S_w(p_n) .* divU(u - u_n);

% Mass conservation (Flow contribution)
pEq2 = @(p, p_n, p_m, tau, sourceFlow)  ...
       flow.alpha_biot .^2 .* S_w(p_n) .* compat(S_w(p_n) .* p_n) ...
       + V .* xi(p_n) .* (p - p_n) ...
       + V .* chi(p_n) .* (S_w(p_m) + C_S(p_m) .* (p - p_m) - S_w(p_n)) ...
       + tau .* divF(Q(p, p_m)) - V .* tau .* sourceFlow;

% Storing equations
model = struct(); % initializing structure to store function handles
model.body = body; % body forces
model.sca2vec = sca2vec; % map from scalar to vector field
model.T = T; % traction forces
model.uEq1 = uEq1; % displacement contribution to momentum equation
model.uEq2 = uEq2; % pressure contribution to momentum equation
model.Q = Q; % darcy flux
model.pEq1 = pEq1; % displacement contribution to mass conservation
model.pEq2 = pEq2; % pressure contribution to mass conservation
model.S_w = S_w; % Water saturation
model.krw = krw; % Relative permeability
model.C_S = C_S; % Specific saturation capacity

end