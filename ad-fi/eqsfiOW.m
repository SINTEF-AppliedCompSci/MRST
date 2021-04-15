function [eqs, state, hst] = eqsfiOW(state0, state, dt, G, W, system, f, varargin)
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
pBH = vertcat(state.wellSol.bhp);
qWs  = vertcat(state.wellSol.qWs);
qOs    = vertcat(state.wellSol.qOs);
%mixs  = vertcat(state.wellSol.mixs);
%mixWs = mixs(:,1);

% previous variables ------------------------------------------------------
p0   = state0.pressure;
sW0  = state0.s(:,1);
pBH0 = vertcat(state0.wellSol.bhp);
qWs0 = vertcat(state0.wellSol.qWs);
qOs0 = vertcat(state0.wellSol.qOs);
%--------------------------------------------------------------------------


%Initialization of independent variables ----------------------------------

if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [p, sW, qWs, qOs, pBH] = ...
            initVariablesADI(p, sW, qWs, qOs, pBH);
    else
        [p0, sW0, tmp, tmp, tmp] = ...
            initVariablesADI(p0, sW0,          ...
            zeros(size(qWs0)), ...
            zeros(size(qOs0)), ...
            zeros(size(pBH0)));                          %#ok
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
bW     = f.bW(p-pcOW);
rhoW   = bW.*f.rhoWS;
% rhoW on face, avarge of neighboring cells (E100, not E300)
rhoWf  = s.faceAvg(rhoW);
%mobW   = trMult.*krW./f.muW(p);
mobW   = trMult.*krW./f.muW(p-pcOW);
dpW     = s.grad(p-pcOW) - g*(rhoWf.*s.grad(G.cells.centroids(:,3)));
% water upstream-index
upc = (value(dpW)>=0);
bWvW = s.faceUpstr(upc, bW.*mobW).*trans.*dpW;


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
    bOvO   = s.faceUpstr(upc, mobO).*trans.*dpO;
else
    mobO   = trMult.*krO./f.muO(p);
    bOvO   = s.faceUpstr(upc, bO.*mobO).*trans.*dpO;
end


% EQUATIONS ---------------------------------------------------------------
% oil:
eqs{1} = (s.pv/dt).*( pvMult.*bO.*(1-sW) - pvMult0.*f.bO(p0).*(1-sW0) ) + s.div(bOvO);


% water:
eqs{2} = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*f.bW(p0).*sW0 ) + s.div(bWvW);


% well equations
if ~isempty(W)
    if ~opt.reverseMode
        wc    = vertcat(W.cells);
        pw   = p(wc);
        rhos = [f.rhoWS, f.rhoOS];
        bw   = {bW(wc), bO(wc)};
        rw   = {};
        mw   = {mobW(wc), mobO(wc)};
        optloc = {'iteration', opt.iteration, ...
                  'model', 'OW', ...
                  'allowWellSignChange', system.well.allowWellSignChange, ...
                  'allowControlSwitching', system.well.allowControlSwitching};
        
        [eqs(3:5), cqs, state.wellSol] = getWellContributions(W, state.wellSol, pBH, {qWs, qOs}, ...
                                                                 pw, rhos, bw, rw, rw, mw, ...
                                                                 optloc{:});

        [wc, cqs] = checkForRepititions(wc, cqs);
        eqs{1}(wc) = eqs{1}(wc) - cqs{2};
        eqs{2}(wc) = eqs{2}(wc) - cqs{1};
    else
        % in reverse mode just gather zero-eqs of correct size
        for eqn = 3:5
            nw = numel(state.wellSol);
            zw = double2ADI(zeros(nw,1), p0);
            eqs(3:5) = {zw, zw, zw};
        end
    end
else % no wells
    eqs(3:5) = {pBH, pBH, pBH};  % empty  ADIs
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





