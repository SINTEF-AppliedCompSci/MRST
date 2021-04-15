function [eqs, state, hst] = eqsfiOWFullT(state0, state, dt, G, W, s, f, varargin)
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
             'minerals', false, ...
             'iteration', -1, ...
             'stepOptions', [],...
             'addfluxes',true);  % Compatibility only

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
T    = state.T;
pBH = vertcat(state.wellSol.bhp);
qWs  = vertcat(state.wellSol.qWs);
qOs    = vertcat(state.wellSol.qOs);
%mixs  = vertcat(state.wellSol.mixs);
%mixWs = mixs(:,1);

% previous variables ------------------------------------------------------
p0  = state0.pressure;
pO0=p0;
sW0 = state0.s(:,1);
T0  = state0.T;
%%
pcOW0 = 0;
if isfield(f, 'pcOW')
    pcOW0  = f.pcOW(sW0);
end
pW0=p0-pcOW0;
assert(all(pW0>0 ))

%--------------------------------------------------------------------------


%Initialization of independent variables ----------------------------------

if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [p, sW, T, qWs, qOs, pBH] = ...
            initVariablesADI(p, sW, T, qWs, qOs, pBH);
    else
        [p0, sW0, T0, tmp, tmp, tmp] = ...
            initVariablesADI(p0, sW0,  T0,        ...
            zeros(size(qWs)), ...
            zeros(size(qOs)), ...
            zeros(size(pBH)));                          %#ok
    end
end
clear tmp
g  = norm(gravity);

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
transMult=1;
if isfield(f, 'transMult')
   transMult=f.transMult(p); 
end
%check for capillary pressure (p_cow)

pcOW = 0;
if isfield(f, 'pcOW')
    pcOW  = f.pcOW(sW);
end

trans=s.T.*transMult;
% -------------------------------------------------------------------------
[krW, krO] = f.relPerm(sW);
%krW = f.krW(sW);
%krO = f.krO(1-sW);

% water props (calculated at oil pressure OK?)
%bW     = f.bW(p);
pW=p-pcOW;
bW     = f.bW(pW,T);
rhoW   = bW.*f.rhoWS;
% rhoW on face, avarge of neighboring cells (E100, not E300)
rhoWf  = s.faceAvg(rhoW);
%mobW   = trMult.*krW./f.muW(p);
mobW   = trMult.*krW./f.muW(pW,T);
dpW     = s.grad(p-pcOW) - g*(rhoWf.*s.grad(G.cells.centroids(:,3)));
% water upstream-index
upc = (double(dpW)>=0);
bWvW = s.faceUpstr(upc, bW.*mobW).*trans.*dpW;


% oil props
bO     = f.bO(p,T);
rhoO   = bO.*f.rhoOS;
rhoOf  = s.faceAvg(rhoO);
dpO    = s.grad(p) - g*(rhoOf.*s.grad(G.cells.centroids(:,3)));
% oil upstream-index
upc = (double(dpO)>=0);

mobO   = trMult.*krO./f.muO(p,T);
bOvO   = s.faceUpstr(upc, bO.*mobO).*trans.*dpO;



% %WELLS ----------------------------------------------------------------
% bWw     = bW(wc);
% bOw     = bO(wc);
% mobWw  = mobW(wc);
% mobOw  = mobO(wc);
%
% %producer mobility
% bWmobWw  = bWw.*mobWw;
% bOmobOw  = bOw.*mobOw;
%
% %set water injector mobility: mobw = mobw+mobo+mobg, mobo = 0;
%
% bWmobWw(iInxW) = bWw(iInxW).*(mobWw(iInxW) + mobOw(iInxW));
% bOmobOw(iInxW) = 0;
% bWmobWw(iInxO) = 0;
% bOmobOw(iInxO) = bOw(iInxO).*(mobWw(iInxO) + mobOw(iInxO));
%
% pw  = p(wc);
%
% bWqW  = -bWmobWw.*Tw.*(pBHP(perf2well) - pw + 0.0*pcOWw + g*dzw.*rhoW(wc));
% bOqO  = -bOmobOw.*Tw.*(pBHP(perf2well) - pw + g*dzw.*rhoO(wc));

% EQUATIONS ---------------------------------------------------------------
% oil:
sO=(1-sW);sO0=(1-sW0);
bO0=f.bO(p0,T0);
eqs{1} = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) + s.div(bOvO);


% water:
bW0=f.bW(p0,T0);
eqs{2} = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + s.div(bWvW);

% temprature
    rhoS={f.rhoOS,f.rhoWS};
    bFsF={bO.*sO,bW.*sW};
    bFsF0={bO0.*sO0,bW0.*sW0};

    bFvF={bOvO,bWvW};
    eF={f.uO(p,T) , f.uW(pW,T)};
    eF0={f.uO(pO0,T0) , f.uW(pW0,T0)};
    hF={f.hO(p,T), f.hW(pW,T)};

    uR=f.uR(T);uR0=f.uR(T0);
    vQ = s.T_r .* s.grad(T);
    %wT=[W.T]';
    %vQqQ=W.WI_r*(wT-T(wc));
    %vQqQ = 0*wT;
    %bFqF={bOqO,bWqW,bGqG};
    % well contributions is taken at the end
    %eqs{3} = (s.pv/dt).*(  (1-pvMult).*uR-(1-pvMult0).*uR0) + s.div( vQ);
    vol=G.cells.volumes;
    eqs{3} = (1./dt).*((vol-pvMult.*s.pv).*uR-(vol-pvMult0.*s.pv).*uR0) + s.div( vQ);
    %eqs{3}=eqs{3};
    %eqs{3}(wc) = eqs{3}(wc)-vQqQ;
    for i=1:numel(eF)       
        eqs{3}  =  eqs{3} + ((s.pv/dt).*( pvMult.*eF{i}.*rhoS{i}.*bFsF{i} - pvMult0.*eF0{i}.*rhoS{i}.*bFsF0{i} )...
                +  s.div( s.faceUpstr(bFvF{i}>0, rhoS{i}.*hF{i}) .* bFvF{i}));         
    end


% well equations
if ~isempty(W)
    if ~opt.reverseMode
        wc    = vertcat(W.cells);
        pw   = p(wc);
        rhos = [f.rhoWS, f.rhoOS];
        bw   = {bW(wc), bO(wc)};
        rw   = {};
        mw   = {mobW(wc), mobO(wc)};
        [eqs(4:6), cqs, state.wellSol, Rw] = getWellContributions(...
            W, state.wellSol, pBH, {qWs, qOs}, pw, rhos, bw, rw, rw, mw, ...
            'iteration', opt.iteration);

        [wc, cqs] = checkForRepititions(wc, cqs);
        eqs{1}(wc) = eqs{1}(wc) - cqs{2};
        eqs{2}(wc) = eqs{2}(wc) - cqs{1};
    else
        % in reverse mode just gather zero-eqs of correct size
        for eqn = 4:6
            nw = numel(state0.wellSol);
            zw = double2ADI(zeros(nw,1), p0);
            eqs(4:6) = {zw, zw, zw};
        end
    end
else % no wells
    eqs(4:6) = {pBH, pBH, pBH};  % empty  ADIs
end
if(~isempty(W))
bFqF={cqs{2},cqs{1}};
hFwp=cell(2,1);

%for i=1:numel(W)
    %{
    hFwp{1}=[hFwp{1};f.rhoOS*repmat(W(i).hO,numel(W(i).cells),1)];
    hFwp{2}=[hFwp{2};f.rhoWS*repmat(W(i).hW,numel(W(i).cells),1)];
    %}
    
    
%end
Tw=[W.T]';
hFwp{2}=f.rhoWS*f.hW(Rw*pBH,Rw*Tw);
hFwp{1}=f.rhoOS*f.hO(Rw*pBH,Rw*Tw);
HqH=cell(3,1);
for i=1:2;
    HqH{i}  = rhoS{i}.*hF{i}(wc).*bFqF{i};
    ind=bFqF{i}>0;
    HqH{i}(ind)  = hFwp{i}(ind).*bFqF{i}(ind);
end
%% add well contributioin
for i=1:2
    eqs{3}(wc) = eqs{3}(wc)  -  HqH{i};
end
end
eqs{3}=eqs{3}/1e6;
if(opt.addfluxes)
    for i=1:numel(bFvF)
        state.bFvF{i}=double(bFvF{i});
    end
end
end
%--------------------------------------------------------------------------



function [wc, cqs] = checkForRepititions(wc, cqs)
[c, ic, ic] = uniqueStable(wc);                                 %#ok<ASGLU>
if numel(c) ~= numel(wc)
    A = sparse(ic, (1:numel(wc))', 1, numel(c), numel(wc));
    wc = c;
    for k=1:numel(cqs)
        cqs{k} = A*cqs{k};
    end
end
end





