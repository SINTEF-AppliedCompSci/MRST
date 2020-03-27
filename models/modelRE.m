function model = modelRE(G, phys_param, krw_faces, mpfa_discr)

% Grid-related quantities
V = G.cell.volumes;

% Discrete mpfa operators
F       = @(x) mpfa_discr.F * x;          % Flux  
boundF  = @(x) mpfa_discr.boundFlux * x;  % Boundary fluxes
divF    = @(x) mpfa_discr.div * x;        % Divergence

% Physical properties
fluid = phys_param.fluid;
vGM   = phys_param.vGM;

% Soil Water Retention Curves (SWRC)
[theta, krw, C_theta] = vanGenuchtenMualemTheta(vGM.alpha, ...
    vGM.theta_s, vGM.theta_r, vGM.n, vGM.m);

% Darcy Flux
Q = @(psi, psi_m) (fluid.rho .* fluid.g ./ fluid.mu) .* krw_faces(psi_m) .* ...
    (F(psi + zetac) + boundF(bcVal));

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
model.psiEQ = psiEq;        % Discrete Mass conservation