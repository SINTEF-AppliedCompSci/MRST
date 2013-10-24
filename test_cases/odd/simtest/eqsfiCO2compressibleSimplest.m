function [eqs, info] = eqsfiCO2compressibleSimplest(state0, state, dt, G, W, s, f, bc, theta, slopedir)

g_cos_t = norm(gravity) * cos(theta);    
g_sin_t = norm(gravity) * sin(theta);    

% current variables
p = state.pressure;  % pressure _measured at top_
h = state.h;         % interface_z - top_z
H = G.cells.H;       % bottom_z - top_z
T = state.temperature;

 %previous variables
p0 = state0.pressure;
h0 = state0.h;

% Initializing adi-variables
[p, h] = initVariablesADI(p, h);

rhoW  = f.water.rho(p, T);     rho0W = f.water.rho(p0, T);
rhoC  = f.CO2.rho(p, T);       rho0C = f.CO2.rho(p0, T);  
betaC = f.CO2.beta(p, T);      beta0C = f.CO2.beta(p0, T);
bderC = f.CO2.bder(p.val, T);  bder0C = f.CO2.bder(p0, T);

poly = VEpolys; % collection of dimensionless VE-correction polynomials.  NB:
                % simplified treatment for water - a whole column is assumed
                % to have the same compressibility, which means we do not
                % establish corretion polynomials for water, only for CO2.

hpi        = @(h)   poly.intAlpha(h, rhoC,  betaC,  bderC,  theta);
hpi0       = @(h)   poly.intAlpha(h, rho0C, beta0C, bder0C, theta);

A          = @(h)   poly.intSqAlpha(h, rhoC,  betaC,  bderC,  theta);
alpha      = @(h)   poly.alpha(h, rhoC,  betaC,  bderC,  theta);

A_some     = @(h, ix) poly.intSqAlpha(h(ix), rhoC(ix), betaC(ix), bderC(ix), theta);
alpha_some = @(h, ix) poly.alpha(h(ix), rhoC(ix), betaC(ix), bderC(ix), theta);

dh = -s.grad(h);  %@ s.grad gives _negative_ gradient, so '-' here.
delta_rho = rhoW - (rhoC .* alpha(h)); % density of CO2 at CO2/water interface

% computing fluxes
dp = - s.grad(p);  % gradient across internal faces (NB: s.grad is
                   % neg. gradient)

dtop = -s.grad(G.cells.z);
pseudo_g = g_cos_t .* dtop; % used in flux computation of CO2 _and_ water

intInx = (prod(double(G.faces.neighbors), 2)~=0); % determine indices of internal faces (hack)
neigh_cells = G.faces.neighbors(intInx,:);
neigh_dists = row_norms(G.cells.centroids(neigh_cells(:,1),:) - ...
                        G.cells.centroids(neigh_cells(:,2),:));
scaling = neigh_dists ./ row_norms(G.faces.normals(intInx,:));
parall_g = - g_sin_t .* G.faces.normals(intInx,:) * slopedir' .* scaling';

CO2_pseudo_g = s.faceAvg(rhoC) .* pseudo_g;
CO2_parall_g = s.faceAvg(rhoC) .* parall_g;
W_parall_g   = s.faceAvg(rhoW) .* parall_g;

upc = logical(dp - CO2_pseudo_g - CO2_parall_g < 0);

% we divide by H, since H is already 'baked in' the computed transmissibilities
% NB: the CO2 flux is a true mass flux (not divided by 'rho')
fluxCO2 = - s.faceUpstr(upc, rhoC .* A(h) ./ f.CO2.mu ./ H) .* s.T .* ...
          (dp - CO2_pseudo_g - CO2_parall_g);

dp_rgh = (s.faceAvg(alpha(h)) .* (dp - CO2_pseudo_g)) - (g_cos_t * s.faceAvg(delta_rho) .* (dh-dtop));
upc = logical(dp_rgh - W_parall_g < 0);

fluxWater = - s.faceUpstr(upc, rhoW .* (H - h) ./ f.water.mu ./ H) .* s.T .* (dp_rgh - W_parall_g);

% Identifying boundary and computing boundary condition
bc_cells = sum(G.faces.neighbors(bc.face, :), 2); % index + 0, or 0 + index
                                                  % = index
pix = find(strcmp(bc.type, 'pressure'));
fix = find(strcmp(bc.type, 'flux'));


% The flux values computed below have signs consistent with whether flow is
% in or out of the cell in question, _not_ according to the established sign
% convention for faces.  In other words, the fluxes (and pressure gradients)
% are positive if they are directed _out_ of the concerned cell, regardless
% of the established orientation of face normals.  The reason for this
% particular treatment is that the boundary contributions are going to be
% added manually, not passed through the 'div' operator, which ensures
% correct orientation for internal faces.
bdist = row_norms(G.cells.centroids(bc_cells,:) - G.faces.centroids(bc.face,:));
bscale = bdist ./ row_norms(G.faces.normals(bc.face,:));

bc_parall_g = - g_sin_t * G.faces.normals(bc.face,:) * slopedir' * bscale;
norm_inv = find(G.faces.neighbors(bc.face,1) == 0);
bc_parall_g(norm_inv) = -bc_parall_g(norm_inv);

% setting boundary conditions for pressure-imposed faces
if pix
    pfaces = bc.face(pix);
    pcells = bc_cells(pix);
    
    bc_dp = bc.value(pix) - p(pcells);
    bc_dp_co2g = bc_dp - (rhoC(pcells) .* bc_parall_g(pix)');
    bc_flux_co2 = - (rhoC(pcells) .* A_some(h, pcells) ./f.CO2.mu ./H(pcells)).* s.T_all(pfaces).* bc_dp_co2g;
    bc_flux_co2(bc_dp_co2g>0) = 0;

    bc_dp_rgh = alpha_some(h, pcells) .* bc_dp; % for simplicity, we consider dh to be zero at boundary...
    bc_dp_rgh = bc_dp_rgh - (rhoW(pcells) .* bc_parall_g(pix)');

    bc_flux_wat = - rhoW(pcells) .* (H(pcells) - h(pcells))./ H(pcells) ./f.water.mu.*s.T_all(pfaces).* bc_dp_rgh;
    
    bc_flux_wat(bc_dp_rgh>0) = -rhoW(pcells(bc_dp_rgh>0))/f.water.mu.* ...
        s.T_all(pfaces(bc_dp_rgh>0)).*bc_dp_rgh(bc_dp_rgh>0);
    
    info.outflowCO2 = bc_flux_co2;
else
    info.outflowCO2 = [];
end


%% Equations


% co2
eqs{1} = (s.pv/dt).*((rhoC.*hpi(h)) - (rho0C.*hpi0(h0))) + s.div(fluxCO2);
for i = numel(W)
    if strcmp(W(i).type, 'Rate')
        cx = W(i).cells;
        eqs{1}(cx) = eqs{1}(cx) - (W(i).val / numel(cx)) ;
    else
        %@ unimplemented
    end
end
if pix
    eqs{1}(pcells)=eqs{1}(pcells) + bc_flux_co2;
end

% water
eqs{2} = (s.pv / dt) .* ((rhoW.*(H-h) - (rho0W.*(H-h0)))) + s.div(fluxWater);

if pix
    eqs{2}(pcells)=eqs{2}(pcells) + bc_flux_wat;
end

end 
%% end of function eqsfiCO2compressibleSimplest

