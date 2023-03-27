function [eqs, info] = eqsfiCO2simplest(state0, state, dt, G, W, s, f, bc, theta, slopedir)

g_cos_t = norm(gravity) * cos(theta);    
g_sin_t = norm(gravity) * sin(theta);    

rhoW = f.water.rho;
rhoC = f.CO2.rho;

% current variables
p = state.pressure;
h = state.h;
H = G.cells.H;

% previous variables
h0 = state0.h;

% Initializing adi-variables
[p, h] = initVariablesADI(p, h);

dh = -s.grad(h);  %@ s.grad gives _negative_ gradient, so '-' here.
drho = rhoW - rhoC;

% computing fluxes
dp = - s.grad(p);  % gradient across internal faces (NB: s.grad is neg. gradient)

dtop = -s.grad(G.cells.z);

pseudo_g = g_cos_t .* dtop;
intInx = (prod(double(G.faces.neighbors), 2)~=0); % determine indices of internal faces (hack)
neigh_cells = G.faces.neighbors(intInx,:);
neigh_dists = row_norms(G.cells.centroids(neigh_cells(:,1),:) - ...
                        G.cells.centroids(neigh_cells(:,2),:));
scaling = neigh_dists ./ row_norms(G.faces.normals(intInx,:));
parall_g = - g_sin_t .* G.faces.normals(intInx,:) * slopedir' .* scaling';

CO2_pseudo_g = rhoC .* pseudo_g;
CO2_parall_g = rhoC .* parall_g;
W_parall_g   = rhoW .* parall_g;

upc = logical(dp - CO2_pseudo_g - CO2_parall_g < 0);
% we divide by H, since H is already 'baked in' the computed transmissibilities
fluxCO2 = - s.faceUpstr(upc, h ./ f.CO2.mu ./ H) .* s.T .* ...
          (dp - CO2_pseudo_g - CO2_parall_g);

dp_rgh = (dp - CO2_pseudo_g) - g_cos_t * drho * (dh-dtop);
upc = logical(dp_rgh - W_parall_g < 0);
fluxWater = - s.faceUpstr(upc, (H - h) ./ f.water.mu ./ H) .* s.T .* (dp_rgh - W_parall_g);
    

% Identifying boundary and computing boundary condition
bc_cells = sum(G.faces.neighbors(bc.face, :), 2); % index + 0, or 0 + index = index
pix = find(strcmp(bc.type, 'pressure'));
fix = find(strcmp(bc.type, 'flux'));
assert(numel(pix) + numel(fix) == numel(bc_cells));

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

if pix
    pfaces = bc.face(pix);
    pcells = bc_cells(pix);
    
    bc_dp = bc.value(pix) - p(pcells);
    bc_dp_co2g = bc_dp - (rhoC * bc_parall_g(pix)');
    bc_flux_co2 = - h(pcells) ./f.CO2.mu ./H(pcells) .* s.T_all(pfaces) .* bc_dp_co2g;
    bc_flux_co2(bc_dp_co2g>0) = 0;

    bc_dp_rgh = bc_dp -  (rhoW * bc_parall_g(pix)');

    bc_flux_wat = -(H(pcells) - h(pcells)) ./ f.water.mu ./H(pcells) .* s.T_all(pfaces) ...
        .*bc_dp_rgh;
    bc_flux_wat(bc_dp_rgh>0) = -1./f.water.mu.*s.T_all(pfaces(bc_dp_rgh>0)).* ...
        bc_dp_rgh(bc_dp_rgh>0);

    info.outflowCO2 = bc_flux_co2;
else
    info.outflowCO2 = [];
end

    
% bc_flux_wat = - (H(bc_cells) - h(bc_cells))./f.muW ./H(bc_cells) .*s.T_all(bc.face).*bc_dp_rgh;
% %bc_flux_wat(bc_dp>0) = - H(bc_cells(bc_dp>0))./f.muW.*s.T_all(bc.face(bc_dp>0)).*bc_dp(bc_dp>0);
% bc_flux_wat(bc_dp_rgh>0) = - 1./f.muW.*s.T_all(bc.face(bc_dp_rgh>0)).*bc_dp(bc_dp_rgh>0);



%% Equations

eqs{1} = (s.pv/dt).*(h - h0) + s.div(fluxCO2);
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
eqs{2} =  (s.pv / dt) .* ((H-h) - (H-h0)) + s.div(fluxWater);

if pix
    eqs{2}(pcells)=eqs{2}(pcells) + bc_flux_wat;
end

end
