function Hw = computeWellHead(w, alpha, rho)
% compute pressure drop from ref depth to connections of corresponding well
% assumes input of class double (not ADI)
if iscell(rho)
    rho  = cell2mat( cellfun(@double, rho, 'UniformOutput', false) );
end

nperf = numel(w.cells);
if ~isfield(w, 'topo')
    topo = [(0:(nperf-1))', (1:nperf)'];
else
    topo = w.topo;
end

rhoMix = sum(alpha.*rho, 2);

dz        = zeros(size(topo,1),1);
dz(1)     = w.dZ(1);
dz(2:end) = w.dZ(topo(2:end,2)) - w.dZ(topo(2:end,1));

g = norm(gravity);

dp        = zeros(size(topo,1),1);
dp(1)     = g*rhoMix(1)*dz(1);
% dp(2:end) = .5*g*rhoMix(topo(2:end,1)).*dz(1:end-1) + .5*g*rhoMix(topo(2:end,2)).*dz(2:end);
if numel(dp)>1
    dp(2:end) = .5*g*rhoMix(topo(2:end,1)).*dz(1:end-1) + .5*g*rhoMix(topo(2:end,2)).*dz(2:end);
end

Hw    = nan(size(topo, 1), 1);
Hw(1) = dp(1);
for k = 2:numel(Hw)
    Hw(topo(k,2)) = Hw(topo(k,1)) + dp(k);
end
end

