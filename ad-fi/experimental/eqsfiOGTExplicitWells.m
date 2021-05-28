function [eqs, state, hst] = eqsfiOGTExplicitWells(state0, state, dt, G, W, s, f, varargin)
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

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'scaling', [],...
             'resOnly', false,...
             'history', []);

opt = merge_options(opt, varargin{:});

if ~isempty(opt.scaling)
    scalFacs = opt.scaling;
else
    scalFacs.rate = 1; scalFacs.pressure = 1;
end

hst = opt.history;

% current variables: ------------------------------------------------------
p    = state.pressure;
T    = state.T;
sG   = state.s(:,2);
if(isfield(state0,'smax'))
    sGmax=max(state.s(:,2),state0.smax(:,2));
else
    sGmax=[];
end
pBHP = vertcat(state.wellSol.bhp);
qGs  = vertcat(state.wellSol.qGs);
qOs  = vertcat(state.wellSol.qOs);
%eQs  = vertcat(state.wellSol.qGs);


% previous variables ------------------------------------------------------
p0  = state0.pressure;
sG0 = state0.s(:,2);
T0  = state0.T;
%--------------------------------------------------------------------------


%Initialization of independent variables ----------------------------------

zw = zeros(size(pBHP));
if ~opt.resOnly,
   % ADI variables needed since we are not only computing residuals.

   if ~opt.reverseMode,

      [p, sG, T, qGs, qOs, pBHP] = ...
         initVariablesADI(p, sG, T, qGs, qOs, pBHP);

   else

      [p0, sG0, zw, zw, zw] = ...
         initVariablesADI(p0, sG0,          ...
                          zeros(size(qGs)), ...
                          zeros(size(qOs)), ...
                          zeros(size(pBHP)));                          %#ok

   end
end

g  = norm(gravity);
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



% -------------------------------------------------------------------------
% [krW, krO] = f.relPerm(sW);
krG = f.krG(sG,'sGmax',sGmax);
krO = f.krOG(1-sG,'sGmax',sGmax);
%check for capillary pressure (p_cow)
pcOG = 0;
if isfield(f, 'pcOG')
    pcOG  = f.pcOG(sG,'sGmax',sGmax);
    pcOGw = pcOG(wc);
end
eG = f.eG(T);
eO = f.eO(T);
eR = f.eR(T);
% water props (calculated at oil pressure OK?)
%bW     = f.bW(p);
bG     = f.bG(p+pcOG);
rhoG   = bG.*f.rhoGS;
% rhoW on face, avarge of neighboring cells (E100, not E300)
rhoGf  = s.faceAvg(rhoG);
%mobW   = trMult.*krW./f.muW(p);
mobG   = trMult.*krG./f.muG(p+pcOG,T);
if(~any(strcmp(G.type,'topSurfaceGrid')))
    dpG     = s.grad(p+pcOG) - g*(rhoGf.*s.grad(G.cells.centroids(:,3)));
else
    dpG     = s.grad(p+pcOG) - g*(rhoGf.*s.grad(G.cells.z));
end
% water upstream-index
upc = (double(dpG)>=0);
bGvG = s.faceUpstr(upc, bG.*mobG).*s.T.*dpG;
eGvG = s.faceUpstr(upc, eG.*mobG).*s.T.*dpG;

% oil props
bO     = f.bO(p);
rhoO   = bO.*f.rhoOS;
rhoOf  = s.faceAvg(rhoO);
mobO   = trMult.*krO./f.muO(p,T);
if(~any(strcmp(G.type,'topSurfaceGrid')))
    dpO    = s.grad(p) - g*(rhoOf.*s.grad(G.cells.centroids(:,3)));
else
    dpO    = s.grad(p) - g*(rhoOf.*s.grad(G.cells.z));
end
% oil upstream-index
upc = (double(dpO)>=0);
bOvO   = s.faceUpstr(upc, bO.*mobO).*s.T.*dpO;
eOvO   = s.faceUpstr(upc, eO.*mobO).*s.T.*dpO;
dT = s.grad(T);

eRvR  = s.TH.*dT;
%eOvO  = eOvO*0.0;
%eGvG  = eGvG*0.0;
%WELLS ----------------------------------------------------------------
bGw     = bG(wc);
eGw     = eG(wc);
bOw     = bO(wc);
eOw     = eO(wc);

mobGw  = mobG(wc);
mobOw  = mobO(wc);

%producer mobility
bGmobGw  = bGw.*mobGw;
bOmobOw  = bOw.*mobOw;
eOmobOw  = eOw.*mobOw;
eGmobGw  = eGw.*mobGw;
%set water injector mobility: mobw = mobw+mobo+mobg, mobo = 0;
WT=[W.T];
if(false)
    bGmobGw(iInxG) = bGw(iInxG).*(mobGw(iInxG) + mobOw(iInxG));
    bOmobOw(iInxG) = 0;
    bGmobGw(iInxO) = 0;
    bOmobOw(iInxO) = bOw(iInxO).*(mobGw(iInxO) + mobOw(iInxO));
else
   bGmobGw(iInxG) = bGw(iInxG)./1e-3;
   bOmobOw(iInxG) = 0;
   bGmobGw(iInxO) = 0;
   bOmobOw(iInxO) = bOw(iInxO)./1e-3;

   if(false)
    eGmobGw(iInxG) = eGw(iInxG)./1e-3;
    eOmobOw(iInxG) = 0;
    eGmobGw(iInxO) = 0;
    eOmobOw(iInxO) = eOw(iInxO)./1e-3;
   else
    eGmobGw(iInxG) = f.eG(WT(iInxG))./1e-3;
    eOmobOw(iInxG) = 0;
    eGmobGw(iInxO) = 0;
    eOmobOw(iInxO) = f.eO(WT(iInxO))./1e-3;
    end
end

pw  = p(wc);

bGqG  = -bGmobGw.*Tw.*(pBHP(perf2well) - pw + 0.0*pcOGw + 0.0*g*dzw.*rhoG(wc));
bOqO  = -bOmobOw.*Tw.*(pBHP(perf2well) - pw + 0.0*g*dzw.*rhoO(wc));

eOqO  = -eOmobOw.*Tw.*(pBHP(perf2well) - pw + 0.0*g*dzw.*rhoO(wc));
eGqG  = -eGmobGw.*Tw.*(pBHP(perf2well) - pw + 0.0*pcOGw + 0.0*g*dzw.*rhoG(wc));

% EQUATIONS ---------------------------------------------------------------


% oil:
eqs{1} = (s.pv/dt).*( pvMult.*bO.*(1-sG) - pvMult0.*f.bO(p0).*(1-sG0) ) + s.div(bOvO);
eqs{1}(wc) = eqs{1}(wc) + bOqO;

% water:
eqs{2} = (s.pv/dt).*( pvMult.*bG.*sG - pvMult0.*f.bG(p0).*sG0 ) + s.div(bGvG);
eqs{2}(wc) = eqs{2}(wc) + bGqG;

eqs{3} = (s.pv/dt).*( pvMult.*(eG.*sG+eO.*(1-sG)+eR) - pvMult0.*(f.eG(T0).*sG0+f.eO(T0).*(1-sG0)+f.eR(T0))) + s.div(eGvG+eOvO+eRvR);
%eqs{3} = (s.pv/dt).*( (pvMult.*(eG.*sG+eO.*(1-sG)+eR)) - pvMult0.*(f.eG(T0).*sG0+f.eO(T0).*(1-sG0)+f.eR(T0))) + s.div(eRvR);
%eqs{3} = (s.pv/dt).*( (pvMult.*( eR)) - pvMult0.*(f.eR(T0))) + s.div(eRvR);
eqs{3}(wc) = eqs{3}(wc) + eGqG + eOqO;
% well equations
eqs{3}=eqs{3}/10.^6;
zeroW = 0*zw;

eqs{4} = Rw'*bGqG + qGs + zeroW;
eqs{5} = Rw'*bOqO + qOs + zeroW;
%eqs{6} = Rw'*(eOqO+eGqG) + eOs + eGs + zeroW;


% Last eq: boundary cond
eqs{6} = handleBC(W, pBHP, [] , qOs, qGs, scalFacs) + zeroW;
end
%--------------------------------------------------------------------------










