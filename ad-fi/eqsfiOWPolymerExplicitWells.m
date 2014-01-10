function [eqs, hst] = eqsfiOWPolymerExplicitWells(state0, state, dt, G, W, s, f, varargin)
% Generate equations for a Oil-Water-Polymer system.

opt = struct('Verbose'    , mrstVerbose, ...
             'reverseMode', false      , ...
             'scaling'    , []         , ...
             'resOnly'    , false      , ...
             'history'    , []         , ...
             'iteration'  , -1         , ...
             'stepOptions', []);  % Ignored.  Framework compatibility only.

opt = merge_options(opt, varargin{:});

if ~isempty(opt.scaling)
    scalFacs = opt.scaling;
else
    scalFacs.rate = 1; scalFacs.pressure = 1;
end

hst = opt.history;

% current variables: ------------------------------------------------------
p     = state.pressure;
sW    = state.s(:,1);
c     = state.c;
cmax0 = state0.cmax;
cmax  = state.cmax;

pBHP = vertcat(state.wellSol.bhp);
qWs  = vertcat(state.wellSol.qWs);
qOs  = vertcat(state.wellSol.qOs);
% poly_num  = vertcat(state.wellSol.poly);
% poly = poly_num;

wciPoly_num = getWellPolymer(W);
wciPoly = wciPoly_num;
% previous variables ------------------------------------------------------
p0  = state0.pressure;
sW0 = state0.s(:,1);
c0  = state0.c;
%--------------------------------------------------------------------------


%Initialization of independent variables ----------------------------------

zw = zeros(size(pBHP));
if ~opt.resOnly,
   % ADI variables needed since we are not only computing residuals.

   if ~opt.reverseMode,

      [p, sW, c, qWs, qOs, wciPoly, pBHP]  = ...
         initVariablesADI(p, sW, c, qWs, qOs, wciPoly_num, pBHP);

   else

      [p0, sW0, c0, zw, zw, zw, zw] = ...
         initVariablesADI(p0, sW0, c0,              ...
                          zeros(size(qWs)),         ...
                          zeros(size(qOs)),         ...
                          zeros(size(wciPoly_num)), ...
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

%Multipliers due to polymer
mixpar = f.mixPar;
cbar   = c/f.cmax;

a = f.muWMult(f.cmax)^(1-mixpar);

b = 1./(1-cbar+cbar./a);

muWMult = b.*f.muWMult(c).^mixpar;

permRed = 1 + ((f.rrf-1)/f.adsMax)*f.ads(max(c, cmax));
muWMult  = muWMult.*permRed;

% polymer injection well:
cw        = c(wc);
cw(iInxW) = wciPoly;
cbarw     = cw/f.cmax;
% muWMultMax = a + (1-a)*cbarw;

%check for capillary pressure (p_cow)
pcOW = 0;
if isfield(f, 'pcOW')
    pcOW  = f.pcOW(sW);
    pcOWw = pcOW(wc);
end

% -------------------------------------------------------------------------
krW = f.krW(sW);
krO = f.krOW(1-sW);
muW = f.muW(p-pcOW);
muWeff = muWMult.*muW;

% water props (calculated at oil pressure OK?)
%bW     = f.bW(p);
bW     = f.bW(p-pcOW);
rhoW   = bW.*f.rhoWS;
% rhoW on face, avarge of neighboring cells (E100, not E300)
rhoWf  = s.faceAvg(rhoW);


mobW   = trMult.*krW./muWeff;
dpW     = s.grad(p-pcOW) - g*(rhoWf.*s.grad(G.cells.centroids(:,3)));
% water upstream-index
upc = (double(dpW)>=0);
bWvW = s.faceUpstr(upc, bW.*mobW).*s.T.*dpW;

% polymer mobility
mobP = (mobW.*c)./(a + (1-a)*cbar);

bWvP = s.faceUpstr(upc, bW.*mobP).*s.T.*dpW;
% bWvP = s.faceUpstr(upc, bW.*mobW).*s.T.*dpW;

% oil props
bO     = f.bO(p);
rhoO   = bO.*f.rhoOS;
rhoOf  = s.faceAvg(rhoO);
mobO   = trMult.*krO./f.BOxmuO(p);
dpO    = s.grad(p) - g*(rhoOf.*s.grad(G.cells.centroids(:,3)));
% oil upstream-index
upc = (double(dpO)>=0);
bOvO   = s.faceUpstr(upc, bO.*mobO).*s.T.*dpO;


%WELLS ----------------------------------------------------------------
bWw     = bW(wc);
bOw     = bO(wc);
mobWw  = mobW(wc);
mobOw  = mobO(wc);

%producer mobility
bWmobWw  = bWw.*mobWw;
bOmobOw  = bOw.*mobOw;

bWmobPw  = (bWmobWw.*cw)./(a + (1-a)*cbarw);

%set water injector mobility: mobw = mobw+mobo+mobg, mobo = 0;
bWmobWw(iInxW) = bWw(iInxW).*(mobWw(iInxW) + mobOw(iInxW));
bOmobOw(iInxW) = 0;

pw  = p(wc);

% Residual equations for source terms
% Transmissibility and pressure differential is common to water and polymer
% in water wells.
tmp = Tw.*(pBHP(perf2well) - pw + pcOWw + g*dzw.*rhoW(wc));
bWqW  = -bWmobWw.*tmp;
% Polymer in water source terms
bWqP  = -bWmobPw.*tmp;
% Oil source terms
bOqO  = -bOmobOw.*Tw.*(pBHP(perf2well) - pw + g*dzw.*rhoO(wc));

% EQUATIONS ---------------------------------------------------------------


% oil:
eqs{1} = (s.pv/dt).*( pvMult.*bO.*(1-sW) - pvMult0.*f.bO(p0).*(1-sW0) ) + s.div(bOvO);
eqs{1}(wc) = eqs{1}(wc) + bOqO;

% water:
eqs{2} = (s.pv/dt).*( pvMult.*bW.*sW     - pvMult0.*f.bW(p0).*sW0     ) + s.div(bWvW);
eqs{2}(wc) = eqs{2}(wc) + bWqW;

% polymer in water:
poro =  s.pv./G.cells.volumes;
eqs{3} =   (s.pv*(1-f.dps)/dt).*(pvMult.*bW.*sW.*c - pvMult0.*f.bW(p0).*sW0.*c0) ...
         + (s.pv/dt).*(f.rhoR*((1-poro)./poro).*(f.ads(max(c, cmax))-f.ads(cmax0))) ... % Adsorption isotherm effect
         + s.div(bWvP);

eqs{3}(wc) = eqs{3}(wc) + bWqP;

% well equations
zeroW = 0*zw;

eqs{4} = Rw'*bWqW + qWs + zeroW;
eqs{5} = Rw'*bOqO + qOs + zeroW;
% Trivial constraint - this is only to get the adjoint partial derivatives
eqs{6} = wciPoly - wciPoly_num + zeroW(iInxW);

% Last eq: boundary cond
eqs{7} = handleBC(W, pBHP, qWs, qOs, [], scalFacs) + zeroW;
end
%--------------------------------------------------------------------------


function wciPoly = getWellPolymer(W)
    if isempty(W)
        wciPoly = [];
        return
    end
    inj   = vertcat(W.sign)==1;
    polInj = cellfun(@(x)~isempty(x), {W(inj).poly});
    polVal = zeros(nnz(inj), 1);
    polVal(polInj) = vertcat(W(inj(polInj)).poly);
    wciPoly = rldecode(polVal, cellfun(@numel, {W(inj).cells}));
end

