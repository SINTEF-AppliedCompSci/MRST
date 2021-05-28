function [eqs, state, hst] = eqsfiOWExplicitWells(state0, state, dt, G, W, system, f, varargin)
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
             'history', [],...
             'temperature', false,...
             'minerals',false, ...
             'iteration', -1, ... % Compatibility only
             'stepOptions', []);  % Compatibility only

opt = merge_options(opt, varargin{:});

if ~isempty(opt.scaling)
    scalFacs = opt.scaling;
else
    scalFacs.rate = 1; scalFacs.pressure = 1;
end

hst = opt.history;
s = system.s;

% current variables: ------------------------------------------------------
p    = state.pressure;
sW   = state.s(:,1);
pBHP = vertcat(state.wellSol.bhp);
qWs  = vertcat(state.wellSol.qWs);
qOs  = vertcat(state.wellSol.qOs);

% previous variables ------------------------------------------------------
p0  = state0.pressure;
sW0 = state0.s(:,1);
%--------------------------------------------------------------------------

if opt.temperature,
   assert (isfield(state, 'T'), ...
          ['State variable must contain valid field ''T'' ', ...
           'when simulating temperature.']);
   T = state.T;
end

if opt.minerals,
   assert (all(isfield(state, {'I', 'M'})), ...
          ['State variable must contain valid fields ''I'' ', ...
           'and ''M'' when simulating minerals.']);

   splitcols = @(a) mat2cell(a, size(a, 1), ones([1, size(a, 2)]));

   I = splitcols(state.I);
   M = splitcols(state.M);
end

%Initialization of independent variables ----------------------------------

zw = zeros(size(pBHP));
if ~opt.resOnly,
   % ADI variables needed since we are not only computing residuals.

   if ~opt.reverseMode,

      if opt.temperature,
         if opt.minerals,

            [p, sW, qWs, qOs, pBHP, T, I{:}, M{:}] = ...
               initVariablesADI(p, sW, qWs, qOs, pBHP, T, I{:}, M{:});

         else

            [p, sW, qWs, qOs, pBHP, T] = ...
               initVariablesADI(p, sW, qWs, qOs, pBHP, T);

         end
      else
         if opt.minerals,

            [p, sW, qWs, qOs, pBHP, I{:}, M{:}] = ...
               initVariablesADI(p, sW, qWs, qOs, pBHP, I{:}, M{:});

         else

            [p, sW, qWs, qOs, pBHP] = ...
               initVariablesADI(p, sW, qWs, qOs, pBHP);

         end
      end

   else

      [p0, sW0, zw, zw, zw] = ...
         initVariablesADI(p0, sW0,          ...
                          zeros(size(qWs)), ...
                          zeros(size(qOs)), ...
                          zeros(size(pBHP)));                          %#ok

   end
end

g  = norm(gravity);
[Tw, dzw, Rw, wc, perf2well, pInx, iInxW, iInxO] = getWellStuff(W);

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
pcOW = 0; pcOWw = 0;
if isfield(f, 'pcOW')
    pcOW  = f.pcOW(sW);
    pcOWw = pcOW(wc);
end

% -------------------------------------------------------------------------
% [krW, krO] = f.relPerm(sW);
krW = f.krW(sW);
krO = f.krOW(1-sW);

% water props (calculated at oil pressure OK?)
bW     = f.bW(p);
% bW     = f.bW(p-pcOW);
rhoW   = bW.*f.rhoWS;
% rhoW on face, avarge of neighboring cells (E100, not E300)
rhoWf  = s.faceAvg(rhoW);
mobW   = trMult.*krW./f.muW(p);
% mobW   = trMult.*krW./f.muW(p-pcOW);
dpW     = s.grad(p-pcOW) - g*(rhoWf.*s.grad(G.cells.centroids(:,3)));
% water upstream-index
upc = (value(dpW)>=0);
bWvW = s.faceUpstr(upc, bW.*mobW).*s.T.*dpW;


% oil props
bO     = f.bO(p);
rhoO   = bO.*f.rhoOS;
rhoOf  = s.faceAvg(rhoO);
dpO    = s.grad(p) - g*(rhoOf.*s.grad(G.cells.centroids(:,3)));
% oil upstream-index
upc = (value(dpO)>=0);
if isfield(f, 'BOxmuO')
    % mob0 is already multplied with b0
    mobO   = trMult.*krO./f.BOxmuO(p);
    bOvO   = s.faceUpstr(upc, mobO).*s.T.*dpO;
else
    if(opt.temperature)
       mobO   = trMult.*krO./f.muO(p,T); 
    else
        mobO   = trMult.*krO./f.muO(p);
    end
    bOvO   = s.faceUpstr(upc, bO.*mobO).*s.T.*dpO;
end


%WELLS ----------------------------------------------------------------
bWw     = bW(wc);
bOw     = bO(wc);
mobWw  = mobW(wc);
mobOw  = mobO(wc);

%producer mobility
bWmobWw  = bWw.*mobWw;
bOmobOw  = bOw.*mobOw;

%set water injector mobility: mobw = mobw+mobo+mobg, mobo = 0;

bWmobWw(iInxW) = bWw(iInxW).*(mobWw(iInxW) + mobOw(iInxW));
bOmobOw(iInxW) = 0;
bWmobWw(iInxO) = 0;
bOmobOw(iInxO) = bOw(iInxO).*(mobWw(iInxO) + mobOw(iInxO));
pw  = p(wc);

bWqW  = -bWmobWw.*Tw.*(pBHP(perf2well) - pw + 0.0*pcOWw + g*dzw.*rhoW(wc));
bOqO  = -bOmobOw.*Tw.*(pBHP(perf2well) - pw + g*dzw.*rhoO(wc));

% EQUATIONS ---------------------------------------------------------------


% oil:
eqs{1} = (s.pv/dt).*( pvMult.*bO.*(1-sW) - pvMult0.*f.bO(p0).*(1-sW0) ) + s.div(bOvO);
eqs{1}(wc) = eqs{1}(wc) + bOqO;

% water:
eqs{2} = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*f.bW(p0).*sW0 ) + s.div(bWvW);
eqs{2}(wc) = eqs{2}(wc) + bWqW;

% well equations
zeroW = 0*zw;

eqs{3} = Rw'*bWqW + qWs + zeroW;
eqs{4} = Rw'*bOqO + qOs + zeroW;

% Last eq: boundary cond
eqs{5} = handleBC(W, pBHP, qWs, qOs, [], scalFacs) + zeroW;

if(opt.temperature)
    %addTemperature(eqs, state0, T, G,  W, dp, mob, fluid, sF,  bF, bF0, bFmobFw, bFqF, pvMult0, pvMult)
    %addTemperature(eqs, state0, T, G,  W, fluid, dp, mob, sF,  bF, bF0, bFmobFw, bFqF, pvMult0, pvMult)
   dp={dpW,dpO};
   mob={mobW,mobO};
   fluid=  {f.eW, f.eO, f.eR};
   sF = {sW, 1-sW};
   bF =  {bW, bO};
   bF0 =  {f.bW(p0), f.bO(p0)};
   %bFmobFw = {bWmobWw, bOmobOw};
   bFqF = {bWqW, bOqO};
   eqs=   addTemperature(eqs, s,dt, state0, T, G,  W, ...
                         fluid,...
                         dp,...
                         mob,...
                         sF,...
                         bF,...
                         bF0,...
                         bFqF,...
                         pvMult0,...
                         pvMult);
end
if(opt.minerals)
   %dp={dpW};
   %mob={mobW};
   %fluid=  f;
   %sF = {sW};
   %bF =  {bW};
   %bF0 =  {f.bW(p0)};
   %bFmobFw = {bWmobWw, bOmobOw};
   %bFqF = {bWqW};
   %eqs = addMinearals()
   eqs = addMinearals(eqs,  I, M, T, state, state0, G,  W, f, f.bW, s, dt, dpW, mobW, bW, sW, sW0, bWqW, pvMult0, pvMult);
end


end
%--------------------------------------------------------------------------










