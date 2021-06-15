function sol = updateConnectionDP(wellmodel, model, sol)
% Explicit update of hydrostatic pressure difference between bottom hole
% and connections based on phase distrubution alonw well-bore.
%
% SYNOPSIS:
%   sol = updateConnectionDP(wellmodel, model, sol)
%
% PARAMETERS:
%   wellmodel   - Simulation well model.
%   model       - Simulation model.
%   sol         - List of current well solution structures
%
% RETURNS:
%   sol         - Well solution structures with updated field 'cdp'
%
% SEE ALSO:
%   WellModel, computeWellContributionsNew.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
   %}
% Explicit update of hydrostatic pressure difference between bottom hole
% and connections.
% input:
% sol : well-solutions with
W = wellmodel.W;
b = wellmodel.bfactors;
rhos = wellmodel.surfaceDensities;
rMax = wellmodel.maxComponents;


nConn       = cellfun(@numel, {W.cells})'; % # connections of each well
perf2well   = rldecode((1:numel(W))', nConn);

% convert to matrices of doubles
toDb  = @(x)cellfun(@double, x, 'UniformOutput', false);
b     = cell2mat(toDb(b));
rMax  = cell2mat(toDb(rMax));
numPh = size(b,2);

% if nargin < 6
%     model = getModel(size(b,2), size(r,2));
% end
% actPh = getActivePhases(model);
%[isActive, actPh] = model.getActivePhases();


for k = 1:numel(sol);
    s  = sol(k);
    w  = W(k);

    if ~isfield(w, 'topo')
        nperf = numel(w.cells);
        w.topo = [(0:(nperf-1))', (1:nperf)'];
    end

    qs = s.cqs; %volumetric in-flux at standard conds
    perfInx = (perf2well == k);
    bk      = b(perfInx,:);
    if ~isempty(rMax)
        rkMax   = rMax(perfInx,:);
    else
        rkMax   = [];
    end

    C = wb2in(w);            % mapping wb-flux to in-flux
    wbqs  = abs(C\qs);       % solve to get well-bore fluxes at surface conds
    wbqst = sum(wbqs, 2);   % total wb-flux at std conds
    % if flux is zero - just use compi
    zi = wbqst == 0;
    if any( zi )
        wbqs(zi,:)  = ones(nnz(zi),1)*w.compi;
        wbqst(zi,:) = sum(wbqs(zi,:), 2);
    end
    % Compute mixture at std conds:
    mixs = wbqs ./ (wbqst*ones(1,numPh));
    % compute volume ratio Vr/Vs
    volRat = compVolRat(mixs, bk, rkMax, model);
    % Mixture density at connection conds (by using static b's)
    rhoMix = (mixs*rhos(:))./volRat;
    % rhoMix is now density between neighboring segments given by
    % topo(:,1)-topo(:,2) computed by using conditions in well-cell
    % topo(:,2). This is probably sufficiently accurate.

    % get dz between segment nodes and bh-node1
    dpt = [0; w.dZ];
    dz  = diff(dpt);
    g   = norm(gravity);
    ddp = g*rhoMix.*dz; % p-diff between connection neighbors
    % well topology assumes we can traverse from top down, but add a loop
    % just in case crazy ordering.
    cdp    = nan(size(ddp));
    cdp(1) = ddp(1);
    its = 0; maxIts = 100;
    while and(any(isnan(cdp)), its<maxIts)
        its = its +1;
        for cnr = 2:numel(cdp)
            cdp(w.topo(cnr,2)) = cdp(w.topo(cnr,1)) + ddp(cnr);
        end
    end
    if its == maxIts
        error(['Problem with topology for well: ', s.name, '. Segments appear not to be connected'])
    end
    sol(k).cdp = cdp;
end
end


function C = wb2in(w)
    conn = w.topo(2:end, :);
    % Number of connections between perforations
    nconn = size(conn, 1);
    % Number of perforations
    nperf = numel(w.cells);
    
    if nconn + 1 ~= nperf
        warning(['Mismatch between connection count (', num2str(nconn+1),...
                ') and perforation count (', num2str(nperf), '). Well model', ...
                'Does not appear to be a tree.']);
    end

    id = (1:nperf)';
    % First identity, then honor topology.
    ii = [id; conn(:, 1)];
    jj = [id; conn(:, 2)];
    vv = [ones(nperf, 1); -ones(nconn, 1)]; 

    C = sparse(ii, jj, vv, nperf, nperf);
end

function volRat = compVolRat(mixs, b, rMax, model)
%
x = mixs;
dg = isprop(model, 'disgas') && model.disgas;
vo = isprop(model, 'vapoil') && model.vapoil;

if dg || vo
    [~, isgas] = model.getVariableField('sg');
    [~, isoil] = model.getVariableField('so');
    
    both = find(isgas | isoil);
    
    g = mixs(:, isgas);
    o = mixs(:, isoil);
    if isa(model, 'ThreePhaseBlackOilModel')
        % Vapoil/disgas
        gor = abs(g./o);
        gor(isnan(gor)) = inf;
        rs = min(rMax(:,1), gor);
        ogr = abs(o./g);
        ogr(isnan(gor)) = inf;
        rv = min(rMax(:,2), ogr);
        d = 1-rs.*rv;
        x(:,isgas) = (x(:,isgas) - rs.*o)./d;
        x(:,isoil) = (x(:,isoil) - rv.*g)./d;
        x(:,both) = x(:,both).*(x(:,both)>0);
    else
        % Only gas dissolution
        x(:,isgas) = x(:,isgas) - rMax(:,1).*o;
        x(:,isgas) = x(:,isgas).*(x(:,isgas)>0);
    end
end

volRat = sum(x./b ,2);
end



