function [eqs, hst, explTrms] = eqsfiBlackOilExplicitWellsOG(state0, state, dt, G, W, s, f, varargin)
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

opt = struct('Verbose',     mrstVerbose,...
             'reverseMode', false,...
             'scaling',     [],...
             'resOnly',     false,...
             'history',     []);
opt = merge_options(opt, varargin{:});

if ~isempty(opt.scaling)
    scalFacs = opt.scaling;
else
    scalFacs.rate = 1; scalFacs.pressure = 1;
end

hst = opt.history;

% current variables: ------------------------------------------------------
p    = state.pressure;
sG   = state.s(:,2);
if(isfield(state0,'smax'))
    sGmax=max(state.s(:,2),state0.smax(:,2));
else
    sGmax=[];
end
rs   = state.rs;
pBHP = vertcat(state.wellSol.bhp);
qOs  = vertcat(state.wellSol.qOs);
qGs  = vertcat(state.wellSol.qGs);

% previous variables ------------------------------------------------------
p0  = state0.pressure;
sG0 = state0.s(:,2);

assert(all(state0.s(:)>=0))
assert(all(state.s(:)>=0))
rs0 = state0.rs;
%--------------------------------------------------------------------------
rsSat  = f.rsSat(p);
isSat = (sG>0 | rs>=rsSat);% | (1 - sW + sG)  == 0;
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

      [p, sG, rs, qOs, qGs, pBHP] = ...
         initVariablesADI(p, sG, rs, qOs, qGs, pBHP);

   else

      [p0, sG0, rs0, zw, zw, zw] = ...
         initVariablesADI(p0, sG0, rs0,     ...
                          zeros(size(qGs)), ...
                          zeros(size(qOs)), ...
                          zeros(size(pBHP)));                          %#ok

   end
end


g  = norm(gravity);
%[Tw, dzw, Rw, wc, perf2well, pInx, iInxW] = getWellStuffOG(W);
[Tw, dzw, Rw, wc, perf2well, pInx, iInxG, iInxO] = getWellStuffOG(W);
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

%check for capillary pressure (p_cog)
pcOG = 0;
if isfield(f, 'pcOG')
    pcOG  = f.pcOG(sG,'sGmax',sGmax);
    pcOGw = pcOG(wc);
end

% -------------------------------------------------------------------------
%[krW, krO, krG] = f.relPerm(sW, sG);
krO = f.krOG(1-sG,'sGmax',sGmax);
krG = f.krG(sG,'sGmax',sGmax);
% water props (calculated at oil pressure OK?)

% oil props
bO     = f.bO(p, rs, isSat);
rhoO   = bO.*(rs*f.rhoGS + f.rhoOS);
rhoOf  = s.faceAvg(rhoO);
muO=f.muO(p,rs,isSat);
mobO   = trMult.*krO./muO;
dpO    = s.grad(p) - g*(rhoOf.*s.grad(G.cells.centroids(:,3)));
% oil upstream-index
upc = (double(dpO)>=0);
bOvO   = s.faceUpstr(upc, bO.*mobO).*s.T.*dpO;
rsbOvO = s.faceUpstr(upc, rs).*bOvO;

% gas props (calculated at oil pressure OK?)
bG     = f.bG(p+pcOG);
%bG     = f.bG(p+pcOG);
rhoG   = bG.*f.rhoGS;
rhoGf  = s.faceAvg(rhoG);
muG = f.muG(p+pcOG);
mobG = trMult.*krG./muG;
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


if ~opt.resOnly
    rhow          = [rhoO.val(wc), rhoG.val(wc)];
else
    rhow          = [rhoO(wc), rhoG(wc)];
end
[Hw, alpha]   = computeWellHeadOG(W, state.wellSol, rhow);

% pressure drawdown
drawdown = -(pBHP(perf2well) + Hw) + p(wc);
%isInj    = double(drawdown) < 0;

% connection mobilities
bOw = bO(wc);
bGw = bG(wc);
rsw = rs(wc);
if(false)
[mobOw, mobGw, crossFlow] = ...
    computeWellMobilitiesOG(W, qOs, qGs, mobO(wc), mobG(wc), ...
                           bOw, bGw, rsw, isInj);

    if any(crossFlow)
    %    fprintf('Crossflow in %2.0d connections\n', nnz(crossFlow));
    end
end
%set water injector mobility: mobw = mobw+mobo+mobg, mobo = 0;
mobGw  = mobG(wc);
mobOw  = mobO(wc);
if(false)
    mobGw(iInxG) = (mobGw(iInxG) + mobOw(iInxG));
    mobOw(iInxG) = 0;
    mobGw(iInxO) = 0;
    mobOw(iInxO) = (mobGw(iInxO) + mobOw(iInxO));
else
    mobGw(iInxG) = 1./1e-3;%muG(wc(iInxG));
    mobOw(iInxG) = 0;
    mobGw(iInxO) = 0;
    mobOw(iInxO) = 1./1e-3;%muO(wc(iInxO));
    %rs(iInxG)=0;
    %rs(iInxO)=0;
end
assert(all(mobOw>=0))
assert(all(mobGw>=0))

if(false)
    qO = (-Tw).*mobOw.*drawdown;
    qG = (-Tw).*mobGw.*drawdown;
else
  qO = (-Tw).*mobOw.*(-(pBHP(perf2well) + 0.0*Hw) + p(wc));
  qG = (-Tw).*mobGw.*(-(pBHP(perf2well) + 0.0*Hw+0.0*pcOGw) + p(wc));
end

bOqO  = bOw.*qO;
bGqG  = bGw.*qG;

%bWqW  = -bWmobWw.*Tw.*(pBHP(perf2well) - pw + pcOWw + g*dzw.*rhoW(wc));
%bOqO  = -bOmobOw.*Tw.*(pBHP(perf2well) - pw + g*dzw.*rhoO(wc));
%bGqG  = -bGmobGw.*Tw.*(pBHP(perf2well) - pw - pcOGw + g*dzw.*rhoG(wc));

% Compute explicit terms
explTrms.wellFlux = [double(qO), double(qG)];
if ~all(sign(Rw'*sum(explTrms.wellFlux,2))==vertcat(W.sign))
%     warning('Some producers are becoming injectors or vice versa !!!!!!')
%     fprintf('Wells changing roles...\n')
end

% EQUATIONS ---------------------------------------------------------------
isSat0 = (double(sG0)>0);
rsSat  = f.rsSat(p);


% oil:
eqs{1} = (s.pv/dt).*( pvMult.*bO.*(1-sG) - pvMult0.*f.bO(p0,rs0,isSat0).*(1-sG0) ) + s.div(bOvO);
eqs{1}(wc) = eqs{1}(wc) - bOqO;
%eqs{1} = addToVals(eqs{1}, wc, bOqO);
%eqs{2} = addToVals(eqs{2}, wc, bWqW);
% gas:
eqs{2} = (s.pv/dt).*...
        ( pvMult.*(bG.*sG + rs.*bO.*(1-sG) ) -...
          pvMult0.*(f.bG(p0).*sG0 + rs0.*f.bO(p0,rs0,isSat0).*(1-sG0) ) )+ ...
          s.div(bGvG + rsbOvO);
%  eqs{2}(wc) = eqs{2}(wc) - bGqG - rsw.*bOqO;
%eqs{2}(wc) = eqs{2}(wc) - bGqG-rsw.*bOqO;

eqs{2}(wc(iInxG)) = eqs{2}(wc(iInxG)) - bGqG(iInxG);
eqs{2}(wc(pInx)) = eqs{2}(wc(pInx)) - bGqG(pInx) - rsw(pInx).*bOqO(pInx);

% closing eqs:
if(isfield(f,'dis_rate'))
   dis_rate=f.dis_rate.*G.cells.volumes;
   %only disolve flowing else hystereis variable should be included
   %dis_rate=dis_rate.*double((rs.val<rsSat.val) & (sG.val>0));
   dis_rate=dis_rate.*double((rs.val<rsSat.val) & (krG.val>0));
   is_sat_loc=(rs.val>=rsSat.val);
   eqs{3} = (s.pv/dt).*...
        ( pvMult.*(rs.*bO.*(1-sG) ) -...
          pvMult0.*(rs0.*f.bO(p0,rs0,isSat0).*(1-sG0) ) )+ ...
          s.div(rsbOvO);
    %eqs{3}(wc) = eqs{2}(wc) - rsw.*bOqO;
    eqs{3}(wc(pInx)) = eqs{3}(wc(pInx)) - rsw(pInx).*bOqO(pInx);
    eqs{3}=eqs{3}-dis_rate;
    eqs{3}(is_sat_loc) = rs(is_sat_loc) - rsSat(is_sat_loc);
else
    eqs{3} = sG + 0*sG0;
    eqs{3}(isSat) = rs(isSat) - rsSat(isSat);
end

% Well equations
% Force wells to be independent ADI variables
zeroW = 0*zw;
eqs{4} = -Rw'*bOqO + qOs + zeroW;
eqs{5} = -Rw'*(bGqG + rsw.*bOqO) + qGs + zeroW;

% Last eq: boundary cond
eqs{6} = handleBC(W, pBHP, [], qOs, qGs, scalFacs) + zeroW;

end

%--------------------------------------------------------------------------

function [Hw, alpha] = computeWellHeadOG(W, wellSol, rho)
inx = 0;
if ~isfield(wellSol, 'flux') || isempty(wellSol(1).flux)
    for wnr = 1:numel(W)
        wellSol(wnr).flux = ones(numel(W(wnr).cells), 1)*W(wnr).compi(:,[2 3])*W(wnr).sign;
    end
end
rhoMix = cell(numel(W),1);
alpha  = cell(numel(W),1);

% the following could be replaced by e.g. a linear system for multilateral
% wells
for wnr = 1:numel(W)
    wbOut = [0 0];
    alpha{wnr}  = zeros(numel(W(wnr).cells), 2);
    rhoMix{wnr} = zeros(numel(W(wnr).cells), 1);
    for cnr = numel(W(wnr).cells):-1:1
        fluxOut = wellSol(wnr).flux(cnr,:);
        if cnr == 1 && W(wnr).sign==1
            wbIn    = sum(wbOut + fluxOut)*W(wnr).compi(:,[2 3]);
        else
            wbIn    = wbOut + fluxOut;
        end
        totIn   = (wbIn    > 0).*wbIn  + ...
                  (wbOut   < 0).*wbOut + ...
                  (fluxOut < 0).*fluxOut;
        if abs(sum(totIn)) < 1e-6/day
            totIn = W(wnr).compi(:,[2 3]);
        end
        if(sum(totIn)>0)
            a = totIn/sum(totIn);
        else
            a=W(wnr).compi(:,[2 3]);
        end
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

function [mO, mG, crossFlow] = computeWellMobilitiesOG(W, qOs, qGs, mO,  mG,  bO, bG, rs, isInj)
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
if(any(iInx))
    injComp = rldecode(vertcat(W(inj).compi), nPerf(inj));
    mtInj   = mO(iInx) + mG(iInx);
    mO(iInx) = injComp(:,2).*mtInj;
    mG(iInx) = injComp(:,3).*mtInj;
end

crossFlow = and(~iInx,isInj);
if any(crossFlow)
    ix = find(crossFlow);
    perf2well = rldecode((1:numel(W))', nPerf);
    qO = qOs(perf2well(ix))./bO(ix);
    qG = (qGs(perf2well(ix))-rs(ix).*qOs(perf2well(ix)))./bG(ix);
    qt =  qO + qG;

    % don't bother with connections having almost zero flow
%     zeroInx = abs(qt.val)<1e-6/day;
%     ix = ix(~zeroInx);
    if ~isempty(ix)
        fprintf('Crossflow detected in %2.0d connections in wells ', numel(ix));
        fprintf('%d ', unique(perf2well(ix)))
        fprintf('\n')
        mt =  mO(ix) + mG(ix);
        mO(ix) = qO.*mt./qt;
        mG(ix) = qG.*mt./qt;

    end
end
assert(all(mO.val>=0))
assert(all(mG.val>=0))
end



