function [eqs, hst, explTrms] = ...
   eqsfiBlackOilExplicitWells(state0, state, dt, G, W, s, f, varargin)
% Generate equations for a Black Oil system (oil, water and gas with gas dissolved in oil).

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

opt = struct('Verbose',     mrstVerbose, ...
             'reverseMode', false,       ...
             'scaling',     [],          ...
             'resOnly',     false,       ...
             'history',     [],          ...
             'stepOptions', []);

opt = merge_options(opt, varargin{:});

if ~isempty(opt.scaling)
    scalFacs = opt.scaling;
else
    scalFacs.rate = 1; scalFacs.pressure = 1;
end

hst = opt.history;

% current variables: ------------------------------------------------------
p    = state.pressure;
sW   = state.s(:,1);
sG   = state.s(:,3);
rs   = state.rs;
pBHP = vertcat(state.wellSol.bhp);
qWs  = vertcat(state.wellSol.qWs);
qOs  = vertcat(state.wellSol.qOs);
qGs  = vertcat(state.wellSol.qGs);

% previous variables ------------------------------------------------------
p0  = state0.pressure;
sW0 = state0.s(:,1);
sG0 = state0.s(:,3);
rs0 = state0.rs;

%--------------------------------------------------------------------------

%isSat = (sG>0) | (1 - sW - sG)  == 0 | rs >= f.rsSat(p)+sqrt(eps);
%isSat = (sG>0) | (1 - sW + sG)  == 0 | rs >= f.rsSat(p);
%isSat = (sG>0) | sW   == 1 | rs >= f.rsSat(p);
isSat = (sG>0) | (1- sW -sG) < sqrt(eps) | rs >= f.rsSat(p);
if isempty(hst)
    hst.numch = zeros(G.cells.num,1);
else
    hst.numch = hst.numch + double(xor(hst.isSat,isSat));
end
hst.isSat = isSat;

%Initialization of independent variables ----------------------------------

zw = zeros(size(pBHP));
if ~opt.resOnly,
   % ADI variables needed since we are not only computing residuals.

   if ~opt.reverseMode,

      [p, sW, sG, rs, qWs, qOs, qGs, pBHP] = ...
         initVariablesADI(p, sW, sG, rs, qWs, qOs, qGs, pBHP);
      zw = 0 * pBHP;

   else

      [p0, sW0, sG0, rs0, zw, zw, zw, zw] = ...
         initVariablesADI(p0, sW0, sG0, rs0, ...
                          zeros(size(qGs)) , ...
                          zeros(size(qWs)) , ...
                          zeros(size(qOs)) , ...
                          zeros(size(pBHP)));                          %#ok

   end
end


g  = norm(gravity);
[Tw, dzw, Rw, wc, perf2well, pInx, iInxW] = getWellStuff(W);

%--------------------
%check for p-dependent tran mult:
trMult = 1;
if isfield(f, 'tranMultR'), trMult = f.tranMultR(p); end

%check for p-dependent porv mult:
pvMult = 1; pvMult0 = 1;
if isfield(f, 'pvMultR')
    pvMult =  f.pvMultR(p);
    pvMult0 = f.pvMultR(p0);
end

%check for capillary pressure (p_cow)
pcOW = 0;
if isfield(f, 'pcOW')
    pcOW  = f.pcOW(sW);
    pcOWw = pcOW(wc);
end
%check for capillary pressure (p_cog)
pcOG = 0;
if isfield(f, 'pcOG')
    pcOG  = f.pcOG(sG);
    pcOGw = pcOG(wc);
end

% -------------------------------------------------------------------------
[krW, krO, krG] = f.relPerm(sW, sG);

% water props (calculated at oil pressure OK?)
bW     = f.bW(p);
%bW     = f.bW(p-pcOW);
rhoW   = bW.*f.rhoWS;
% rhoW on face, avarge of neighboring cells (E100, not E300)
rhoWf  = s.faceAvg(rhoW);
mobW   = trMult.*krW./f.muW(p);
%mobW   = trMult.*krW./f.muW(p-pcOW);
dpW     = s.grad(p-pcOW) - g*(rhoWf.*s.grad(G.cells.centroids(:,3)));
% water upstream-index
upc = (double(dpW)>=0);
bWvW = s.faceUpstr(upc, bW.*mobW).*s.T.*dpW;


% oil props
bO     = f.bO(p, rs, isSat);
rhoO   = bO.*(rs*f.rhoGS + f.rhoOS);
rhoOf  = s.faceAvg(rhoO);
mobO   = trMult.*krO./f.muO(p,rs,isSat);
dpO    = s.grad(p) - g*(rhoOf.*s.grad(G.cells.centroids(:,3)));
% oil upstream-index
upc = (double(dpO)>=0);
bOvO   = s.faceUpstr(upc, bO.*mobO).*s.T.*dpO;
rsbOvO = s.faceUpstr(upc, rs).*bOvO;

% gas props (calculated at oil pressure OK?)
bG     = f.bG(p);
%bG     = f.bG(p+pcOG);
rhoG   = bG.*f.rhoGS;
rhoGf  = s.faceAvg(rhoG);
mobG = trMult.*krG./f.muG(p);
%mobG = trMult.*krG./f.muG(p+pcOG);

dpG     = s.grad(p+pcOG) - g*(rhoGf.*s.grad(G.cells.centroids(:,3)));
% water upstream-index
upc = (double(dpG)>=0);
bGvG = s.faceUpstr(upc,bG.*mobG).*s.T.*dpG;

%WELLS ----------------------------------------------------------------
%mobWw  = mobW(wc);
%mobOw  = mobO(wc);
%mobGw  = mobG(wc);

% well bore head term computed explicitly, also compute connection phase
% volume fractions (alpha)


if ~opt.resOnly && ~opt.reverseMode
    rhow          = [rhoW.val(wc), rhoO.val(wc), rhoG.val(wc)];
else
    rhow          = [rhoW(wc), rhoO(wc), rhoG(wc)];
end
if(~isempty(W))
    [Hw, alpha]   = computeWellHead(W, state.wellSol, rhow);
else
    Hw=0*bW(wc);
end
% pressure drawdown
drawdown = -(pBHP(perf2well) + Hw) + p(wc);
isInj    = double(drawdown) < 0;

% connection mobilities
bWw = bW(wc);
bOw = bO(wc);
bGw = bG(wc);
rsw = rs(wc);
if(~isempty(W))
[mobWw, mobOw, mobGw, crossFlow] = ...
    computeWellMobilities(W, qWs, qOs, qGs, mobW(wc), mobO(wc), mobG(wc), ...
                          bWw, bOw, bGw, rsw, isInj);

if any(crossFlow)
%    fprintf('Crossflow in %2.0d connections\n', nnz(crossFlow));
end
else
    mobWw=0*bOw;
    mobOw=0*bOw;
    mobGw=0*bOw;
end
qW = (-Tw).*mobWw.*drawdown;
qO = (-Tw).*mobOw.*drawdown;
qG = (-Tw).*mobGw.*drawdown;


bWqW  = bWw.*qW;
bOqO  = bOw.*qO;
bGqG  = bGw.*qG;

%bWqW  = -bWmobWw.*Tw.*(pBHP(perf2well) - pw + pcOWw + g*dzw.*rhoW(wc));
%bOqO  = -bOmobOw.*Tw.*(pBHP(perf2well) - pw + g*dzw.*rhoO(wc));
%bGqG  = -bGmobGw.*Tw.*(pBHP(perf2well) - pw - pcOGw + g*dzw.*rhoG(wc));

% Compute explicit terms
explTrms.wellFlux = [double(qW), double(qO), double(qG)];
if(~isempty(W))
if ~all(sign(Rw'*sum(explTrms.wellFlux,2))==vertcat(W.sign))
%     warning('Some producers are becoming injectors or vice versa !!!!!!')
%     fprintf('Wells changing roles...\n')
end
end

% EQUATIONS ---------------------------------------------------------------
isSat0 = (double(sG0)>0);
rsSat  = f.rsSat(p);


% oil:
eqs{1} = (s.pv/dt).*( pvMult.*bO.*(1-sW-sG) - pvMult0.*f.bO(p0,rs0,isSat0).*(1-sW0-sG0) ) + s.div(bOvO);
eqs{1}(wc) = eqs{1}(wc) - bOqO;
%eqs{1} = addToVals(eqs{1}, wc, bOqO);
% water:
eqs{2} = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*f.bW(p0).*sW0 ) + s.div(bWvW);
eqs{2}(wc) = eqs{2}(wc) - bWqW;
%eqs{2} = addToVals(eqs{2}, wc, bWqW);
% gas:
eqs{3} = (s.pv/dt).*...
        ( pvMult.*(bG.*sG + rs.*bO.*(1-sW-sG) ) -...
          pvMult0.*(f.bG(p0).*sG0 + rs0.*f.bO(p0,rs0,isSat0).*(1-sW0-sG0) ) )+ ...
          s.div(bGvG + rsbOvO);
eqs{3}(wc) = eqs{3}(wc) - bGqG - rsw.*bOqO;

% closing eqs:
eqs{4} = sG + 0*sG0;
eqs{4}(isSat) = rs(isSat) - rsSat(isSat);

% Well equations
% Force wells to be independent ADI variables
zeroW = 0*zw;

eqs{5} = -Rw'*bWqW + qWs + zeroW;
eqs{6} = -Rw'*bOqO + qOs + zeroW;
eqs{7} = -Rw'*(bGqG + rsw.*bOqO) + qGs + zeroW;

% Last eq: boundary cond
if(~isempty(W))
    eqs{8} = handleBC(W, pBHP, qWs, qOs, qGs, scalFacs) + zeroW;
else
    eqs{8} = zeroW;
end
for i=1:numel(eqs)
   assert(all(isfinite(double(eqs{i}))));
end

end

%--------------------------------------------------------------------------

function [Hw, alpha] = computeWellHead(W, wellSol, rho)
inx = 0;
for wnr = 1:numel(W)
    if ~isfield(wellSol(wnr), 'flux') || isempty(wellSol(wnr).flux)
        wellSol(wnr).flux = ones(numel(W(wnr).cells), 1)*W(wnr).compi*W(wnr).sign;
    end
end

rhoMix = cell(numel(W),1);
alpha  = cell(numel(W),1);

% the following could be replaced by e.g. a linear system for multilateral
% wells
for wnr = 1:numel(W)
    wbOut = [0 0 0];
    alpha{wnr}  = zeros(numel(W(wnr).cells), 3);
    rhoMix{wnr} = zeros(numel(W(wnr).cells), 1);
    for cnr = numel(W(wnr).cells):-1:1
        fluxOut = wellSol(wnr).flux(cnr,:);
        if cnr == 1 && W(wnr).sign==1
            wbIn    = sum(wbOut + fluxOut)*W(wnr).compi;
        else
            wbIn    = wbOut + fluxOut;
        end
        totIn   = (wbIn    > 0).*wbIn  + ...
                  (wbOut   < 0).*wbOut + ...
                  (fluxOut < 0).*fluxOut;
        if abs(sum(totIn)) < 1e-6/day
            totIn = W(wnr).compi;
        end
        a = totIn/sum(totIn);
        alpha{wnr}(cnr,:)  = a;
        rhoMix{wnr}(cnr) = sum(a.*rho(inx+cnr,:));
        wbOut = wbIn;
    end
    inx = inx + numel(W(wnr).cells);
end
rhoMixAvg = cellfun(@(x)[x(1); .5*(x(1:end-1)+x(2:end))], rhoMix, 'UniformOutput', false);

g   = norm(gravity);
dzw = arrayfun(@(x)[x.dZ(1); x.dZ(2:end)-x.dZ(1:end-1)], W, 'UniformOutput', false);
ddp = cellfun(@(x,y)g*x.*y, dzw, rhoMixAvg, 'UniformOutput', false);
Hw  = cellfun(@cumsum, ddp, 'UniformOutput', false);
Hw = vertcat(Hw{:});
alpha = vertcat(alpha{:});
%dzw = vertcat(W.dZ);
%dpw = g*dzw.*vertcat(rhoMixAvg{:});
end

%--------------------------------------------------------------------------

function [mW, mO, mG, crossFlow] = computeWellMobilities(W, qWs, qOs, qGs, mW, mO, mG, bW, bO, bG, rs, isInj)
% for producer producing connections, m = [mW mO mG] remains unchanged
% for producer injecting connections the (assumed uniform) flow in the
% well-bore is calculated by
%   qW = qWs/bW
%   qO = qOs/bO
%   qG = (qGs-rs*qOs)/bG
% and the mobility as m = [qW qO qG]*(mW+mO+mG)/(qW+qO+qG)
%
% for injector injecting connections the mobility is calculated by
% m = compi*(mW+mO+mG)
% injector producing connections is treated as above, but supresses a
% warning

nPerf = arrayfun(@(x)numel(x.cells), W)';
inj   = vertcat(W.sign)==1;
iInx  = rldecode(inj, nPerf);
% update crossflow flag
%check for injector producing connections
% if any(and(iInx, ~isInj))
%     warning('Crossflow detected in injectors, no special treatment for this case')
% end
injComp = rldecode(vertcat(W(inj).compi), nPerf(inj));
if ~isempty(injComp)
mtInj   = mW(iInx) + mO(iInx) + mG(iInx);
mW(iInx) = injComp(:,1).*mtInj;
mO(iInx) = injComp(:,2).*mtInj;
mG(iInx) = injComp(:,3).*mtInj;
end
crossFlow = and(~iInx,isInj);
if any(crossFlow)
    ix = find(crossFlow);
    perf2well = rldecode((1:numel(W))', nPerf);

    qW = qWs(perf2well(ix))./bW(ix);
    qO = qOs(perf2well(ix))./bO(ix);
    qG = (qGs(perf2well(ix))-rs(ix).*qOs(perf2well(ix)))./bG(ix);
    qt = qW + qO + qG;

    % don't bother with connections having almost zero flow
%     zeroInx = abs(qt.val)<1e-6/day;
%     ix = ix(~zeroInx);
    if ~isempty(ix)
        fprintf('Crossflow detected in %2.0d connections in wells ', numel(ix));
        fprintf('%d ', unique(perf2well(ix)))
        fprintf('\n')
        mt = mW(ix) + mO(ix) + mG(ix);
        mW(ix) = qW.*mt./qt;
        mO(ix) = qO.*mt./qt;
        mG(ix) = qG.*mt./qt;

    end
end
end



