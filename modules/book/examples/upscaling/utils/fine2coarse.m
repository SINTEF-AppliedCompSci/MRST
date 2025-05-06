function sol = fine2coarse(solf, G, W)
% project fine scale solution to coarse grid by volume averageing pressure
% and summing fluxes

% coarse grid pressure
volf = G.parent.cells.volumes;
pc = accumarray(G.partition, solf.pressure.*volf)./...
     accumarray(G.partition, volf); % = Gc.volumes

np = size(solf.s,2);
sc = zeros(G.cells.num, np);
for k = 1:np
    sc(:,k) = accumarray(G.partition, solf.s(:,k).*volf)./...
              accumarray(G.partition, volf); % = Gc.volumes
end
% coarse grid fluxes
fluxc = coarsenFlux(G, solf.flux);

% coarse grid face pressure:
% Since we have no half-trans, there is no unambigous choice for interior
% face pressure, pressure for exterior faces can be choosen 'freely', e.g.,
% use areal averaging
if isfield(solf, 'facePressure'),
    sfp  = solf.facePressure(G.faces.fconn); % sub face pressures
    sfa  = G.parent.faces.areas(G.faces.fconn);   % sub face areas
    % standard incantation:
    fcno = rldecode(1:G.faces.num, diff(G.faces.connPos), 2)';
    % finally coarse face pressures:
    pfc  = accumarray(fcno, sfp.*sfa)./accumarray(fcno, sfa);
    % set interior face pressures to nan to avoid using these for anything
    pfc( all(G.faces.neighbors, 2) ) = nan;
else
    pfc = nan(G.faces.num, 1);
end

% coarse well fluxes (bhp remain the same)
if ~isempty(W),
    wsc = solf.wellSol;
    for i = 1 : numel(W)
        % standard ...
        pno = rldecode(1 : numel(W(i).cells), diff(W(i).fcellspos), 2).';
        wsc(i).flux = accumarray(pno, solf.wellSol(i).flux(W(i).fperf));
    end
else
    wsc = [];
end

% update solc
sol = struct('pressure',      pc,     ...
              'flux',          fluxc,         ...
              'facePressure',  pfc, ...
              's',              sc, ...
              'wellSol',       wsc);
end
