function model = modelRE(G, phys, krw_faces, mpfa_discr, bcVal, gEffects)

% Gravity effects
if strcmp(gEffects, 'on')
    gravOn = 1;
else
    gravOn = 0;
end

% Grid-related quantities
zc = G.cells.centroids(:, end);  % cell centers in z-direction
zf = G.faces.centroids(:, end);  % face centers in z-direction
zetac = max(zf) - zc;            % centroids of cells of elev. head
V = G.cells.volumes;             % Cell volumes

% Discrete mpfa operators
F       = @(x) mpfa_discr.F * x;          % Flux  
boundF  = @(x) mpfa_discr.boundFlux * x;  % Boundary fluxes
divF    = @(x) mpfa_discr.div * x;        % Divergence

% Soil Water Retention Curves (SWRC)
[theta, krw, C_theta] = vanGenuchtenMualemTheta(phys.flow.alpha, ...
    phys.flow.theta_s, phys.flow.theta_r, phys.flow.n, phys.flow.m);

% Darcy Flux
Q = @(psi, psi_m) (phys.flow.gamma ./ phys.flow.mu) .* krw_faces(psi_m) .* ...
    (F(psi + gravOn * zetac) + boundF(bcVal));

% Mass Conservation                         
psiEq = @(psi, psi_n, psi_m, tau, source)   (V./tau) .* (theta(psi_m) ...
    + C_theta(psi_m) .* (psi - psi_m) - theta(psi_n)) ...
    + divF(Q(psi, psi_m)) - V .* source;

% Model equations
model = [];
model.theta = theta;        % Consitutive relationship for water content
model.krw = krw;            % Consitutive relationship for rel perm
model.C_theta = C_theta;    % Consitutive relationship specific moisture capacity
model.Q = Q;                % Discrete Darcy equation
model.psiEq = psiEq;        % Discrete Mass conservation