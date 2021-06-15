function sol = updateConnDP(W, sol, b, rMax, rhos, model)
% Explicit update of hydrostatic pressure difference between bottom hole
% and connections.
% input:
% sol : well-solutions with

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

nConn       = cellfun(@numel, {W.cells})'; % # connections of each well
perf2well   = rldecode((1:numel(W))', nConn);

% convert to matrices of doubles
toDb  = @(x)cellfun(@value, x, 'UniformOutput', false);
b     = cell2mat(toDb(b));
rMax  = cell2mat(toDb(rMax));
numPh = size(b,2);

if nargin < 6
    model = getModel(size(b,2), size(r,2));
end
actPh = getActivePhases(model);

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
        wbqs(zi,:)  = ones(nnz(zi),1)*w.compi(actPh);
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

function model = getModel(state)
np = size(state.s, 2); % #states
nr = 0; % #solutions
if isfield(state, 'rs'), nr = nr+1;end
if isfield(state, 'rv'), nr = nr+1;end
models = {''  , ''  , ''
          'OW', ''  , ''
          '3P', 'BO', 'VO'};
model = models{np, nr+1};
end


% function qs = surfrates(bq, r, model)
% qs = bq;
% switch model
%     case 'OW' %done
%     case '3P' %done
%     case 'BO'
%         qs(:,3) = qs(:,3) + r(:,1)*bq(:,2);
%     case 'VO'
%         qs(:,3) = qs(:,3) + r(:,1)*bq(:,2);
%         qs(:,2) = qs(:,2) + r(:,2)*bq(:,3);
%     otherwise
%         error(['Unknown model: ', model]);
% end
% end

% function vc = surf2conn(vs, b, r, model)
% x = vs;
% switch model
%     case 'OW' %done
%     case '3P' %done
%     case 'BO'
%         x(:,3) = x(:,3) - r(:,1)*vs(:,2);
%     case 'VO'
%         x(:,3) = x(:,3) - r(:,1)*vs(:,2);
%         x(:,2) = x(:,2) - r(:,2)*vs(:,3);
%     otherwise
%         error(['Unknown model: ', model]);
% end
% vc = x./b;
% end

function C = wb2in(w)
    nperf = numel(w.cells);
    ii = [w.topo(:,2); w.topo(2:end, 1)];
    jj = [(1:nperf)'; (2:nperf)'];
    vv = [ones(nperf, 1); -ones(nperf-1, 1)];
    C = sparse(ii, jj, vv, nperf, nperf);
end

function volRat = compVolRat(mixs, b, rMax, model)
%
x = mixs;
switch model
    case 'OW' %done
    case 'WG' %done
    case '3P' %done
    case 'BO'
        x(:,3) = x(:,3) - rMax(:,1).*mixs(:,2);
        x(:,3) = x(:,3).*(x(:,3)>0);
    case 'VO'
        gor = abs(mixs(:,3)./mixs(:,2));
        gor(isnan(gor)) = inf;
        rs = min(rMax(:,1), gor);
        ogr = abs(mixs(:,2)./mixs(:,3));
        ogr(isnan(gor)) = inf;
        rv = min(rMax(:,2), ogr);
        d = 1-rs.*rv;
        x(:,3) = (x(:,3) - rs.*mixs(:,2))./d;
        x(:,2) = (x(:,2) - rv.*mixs(:,3))./d;
        x(:,2:3) = x(:,2:3).*(x(:,2:3)>0);
    otherwise
        error(['Unknown model: ', model]);
end
volRat = sum(x./b ,2);
end



