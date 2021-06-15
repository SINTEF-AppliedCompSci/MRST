function [eqs, bWqW, bOqO, bGqG, sol] = eqsWellBO(w, pBH, qs, p, rho, b, rs, m, sol)
%Undocumented Utility Function

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

if ~isfield(sol, 'alpha')
    alpha = [0 0 0];
else
    alpha = sol.alpha;
end
nperf = numel(w.cells);
if ~isfield(w, 'topo')
    topo = [(0:(nperf-1))', (1:nperf)'];
else
    topo = w.topo;
end
nperf = numel(w.cells);

qWs = qs{1}; qOs = qs{2}; qGs = qs{3};
bW  = b{1};  bO  = b{2};  bG  = b{3};
mW  = m{1};  mO  = m{2};  mG  = m{3};

Hw       = computeWellHead(w, alpha, rho);
seg_pres = double(pBH)+Hw;
drawdown = -pBH - Hw + p;

if sol.sign < 0 % producer
    crossflow = double(drawdown)<0;
    if any(crossflow)
        inj = find(crossflow);
        ni  = numel(inj);
        % compute total production at injection connection conditions
        qW = qWs./bW(inj);
        qO = qOs./bO(inj);
        qG = (qGs - rs(inj)*qOs)./bG(inj);
        qt = qW + qO + qG;

        mt = mW(inj) + mO(inj) + mG(inj);
        mW(inj) = qW.*mt./qt;
        mO(inj) = qO.*mt./qt;
        mG(inj) = qG.*mt./qt;
    end
elseif sol.sign > 0 % injector
    crossflow = double(drawdown)>0;
    if ~any(crossflow)
        mt  = mW + mO + mG;
        mW  = w.compi(1)*mt;
        mO  = w.compi(2)*mt;
        mG  = w.compi(3)*mt;
    else
        prd = find(crossflow);
        %compute total in-flux into wellbore at surf conditions
        bWqW = -w.WI(prd).*bW(prd).*mW(prd).*drawdown(prd);
        bOqO = -w.WI(prd).*bO(prd).*mO(prd).*drawdown(prd);
        bGqG = -w.WI(prd).*bG(prd).*mG(prd).*drawdown(prd);

        qInWs = qWs - sum(bWqW);
        qInOs = qOs - sum(bOqO);
        qInGs = qGs - sum(bGqG + rs(prd).*bOqO);

        %compute total in-flux into wellbore at injection connection conditions
        inj = find(~crossflow);
        qW = qInWs./bW(inj);
        qO = qInOs./bO(inj);
        qG = (qInGs - rs(inj).*qInOs)./bG(inj);
        qt = qW + qO + qG;

        mt = mW(inj) + mO(inj) + mG(inj);
        mW(inj) = qW.*mt./qt;
        mO(inj) = qO.*mt./qt;
        mG(inj) = qG.*mt./qt;
    end
end

bWqW = -w.WI.*bW.*mW.*drawdown;
bOqO = -w.WI.*bO.*mO.*drawdown;
bGqG = -w.WI.*bG.*mG.*drawdown;

eqs{1} = -sum(bWqW) + qWs;
eqs{2} = -sum(bOqO) + qOs;
eqs{3} = -sum(bGqG + rs.*bOqO) + qGs;
eqs{4} = handleBC(sol, pBH, qWs, qOs, qGs);

if nargout > 4
    % finally explicit computations of alpha:
    % 'divergence'-matrix for wellbore:
    ii = [topo(:,2); topo(2:end, 1)];
    jj = [(1:nperf)'; (2:nperf)'];
    vv = [ones(nperf, 1); -ones(nperf-1, 1)];
    C = sparse(ii, jj, vv, nperf, nperf);
    qs   = [double(bWqW), double(bOqO), double(bGqG) + double(rs).*double(bOqO)];
    wbqs = C\qs;
    % take flow in each segment as equal to upstream:
    % reindex to avoid 0-index (producers)
    segInflows = zeros(nperf,3);
    dwnFlow   = sum(wbqs,2)>0;
    segInflows(topo(dwnFlow,2),:) = wbqs(dwnFlow,:);
    upFlow    = [false; sum(wbqs(2:end,:),2)<0];
    segInflows(topo(upFlow,1),:) = -wbqs(upFlow,:);
    isPrd  = sum(qs,2)<0;
    segInflows(isPrd, :) = segInflows(isPrd, :) - qs(isPrd,:);

    % at connection conditions:
    segq   = [segInflows(:,1)./double(bW), ...
        segInflows(:,2)./double(bO), ...
        (segInflows(:,3) - double(rs).*segInflows(:,2))./double(bG)];
    sol.alpha  = segq./(sum(segq,2)*[1 1 1]);

    sol.seg_pres = seg_pres;
end
end


