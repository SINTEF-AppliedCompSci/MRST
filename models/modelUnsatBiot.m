function [model] = modelUnsatBiot(G, phys, krw_faces, mpfa_discr, mpsa_discr, ...
    bcFlowVals, bcMechVals)

% Physical properties
flow = phys.flow;
mech = phys.mech;

% Soil Water Retention Curves (SWRC)
[S_w, krw, C_S] = vanGenuchtenMualemSw(flow.a, flow.S_r, flow.n, flow.m);

% Grid-related quantities
Nc = G.cells.num;     % total number of cells
Nd = G.griddim;       % grid dimension
V  = G.cells.volumes; % cell volumes
zc = G.cells.centroids(:,end); % cell centers in z-direction     
zf = G.faces.centroids(:,end); % face centers in z-direction
zetac = max(zf) - zc; % cell centers of elev. head

% MPFA/MPSA discrete operators

% MPSA operators
S           = @(x) mpsa_discr.stress * x; % stress
boundS      = @(x) mpsa_discr.boundStress * x; % stress boundary
gradP       = @(x) mpsa_discr.gradP * x; % gradient of pressure
divU        = @(x) mpsa_discr.divD * x; % divergence of displacement
compat      = @(x) mpsa_discr.stabDelta * x; % stability operator
divS        = @(x) mpsa_discr.div * x; % divergence of stress

% MPFA operators
F           = @(x) mpfa_discr.F * x; % flux
boundF      = @(x) mpfa_discr.boundFlux * x; % boundary fluxes
divF        = @(x) mpfa_discr.div * x; % divergence of flux

% Mechanical discrete equations

% Function that maps "c" to "Nd*c"
sca2vec = @(x) repmat(x, [Nd, 1]); 

% Body forces
g_body = zeros(Nd * Nc, 1);  % intializing body forces vector
g_body(Nd:Nd:end) = flow.g;  % assigning g to z-cells
body = @(p_n) ((1 - flow.poro) .* mech.rho + ...
    flow.poro .* sca2vec(S_w(p_n)) .* flow.rho) .* g_body;

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
Q = @(p, p_m) (krw_faces(p_m) ./ flow.mu) .* (F(p + flow.gamma .* zetac) ...
    + boundF(bcFlowVals)); 

% Mass conservation (Mechanics contribution)
pEq1 = @(p_n, u, u_n) flow.alpha_biot .* S_w(p_n) .* divU(u - u_n);

% Mass conservation (Flow contribution)
pEq2 = @(p, p_n, p_m, tau, sourceFlow)  ...
       flow.alpha_biot .^2 .* S_w(p_n) .* compat(S_w(p_n) .* p_n) ...
       + V .* xi(p_n) .* (p - p_n) ...
       + V .* chi(p_n) .* (S_w(p_m) + C_S(p_m) .* (p - p_m) - S_w(p_n)) ...
       + tau .* divF(Q(p, p_m)) - V .* sourceFlow;

% Storing equations
model = struct();   % initializing structure
model.body = body;  % body forces
model.T = T;        % traction forces
model.uEq1 = uEq1;  % displacement contribution to momentum equation
model.uEq2 = uEq2;  % pressure contribution to momentum equation
model.Q = Q;        % darcy flux
model.pEq1 = pEq1;  % displacement contribution to mass conservation
model.pEq2 = pEq2;  % pressure contribution to mass conservation
model.S_w = S_w;    % Water saturation
model.krw = krw;    % Relative permeability
model.C_S = C_S;    % Specific saturation capacity

end