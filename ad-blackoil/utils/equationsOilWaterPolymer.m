function [problem, state] = equationsOilWaterPolymer(state0, state, ...
   model, dt, drivingForces, varargin)
% Get linearized problem for oil/water/polymer system with black oil-style
% properties
opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'iteration', -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.Wells;
assert(isempty(drivingForces.bc) && isempty(drivingForces.src))

% Operators, grid and fluid model.
s = model.operators;
G = model.G;
f = model.fluid;

% Properties at current timestep
[p, sW, c, cmax, wellSol] = model.getProps(state, 'pressure', 'water', ...
   'polymer', 'polymermax', 'wellsol');

% Properties at previous timestep
[p0, sW0, c0, cmax0] = model.getProps(state0, 'pressure', 'water', ...
   'polymer', 'polymermax');

pBH    = vertcat(wellSol.bhp);
qWs    = vertcat(wellSol.qWs);
qOs    = vertcat(wellSol.qOs);
qWPoly = vertcat(wellSol.qWPoly);

% Initialize independent variables.
if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [p, sW, c, qWs, qOs, qWPoly, pBH] = ...
            initVariablesADI(p, sW, c, qWs, qOs, qWPoly, pBH);
    else
        [p0, sW0, c0, tmp, tmp, tmp, tmp] = ...
            initVariablesADI(p0, sW0, c0,       ...
            zeros(size(qWs)), ...
            zeros(size(qOs)), ...
            zeros(size(qWPoly)), ...
            zeros(size(pBH)));                          %#ok
        clear tmp
    end
end

% We will solve for pressure, water saturation (oil saturation follows via
% the definition of saturations), polymer concentration and well rates +
% bhp.
primaryVars = {'pressure', 'sW', 'polymer', 'qWs', 'qOs', 'qWPoly', 'bhp'};

% Pressure dependent transmissibility multiplier 
[trMult, pvMult, pvMult0, transMult] = deal(1);
if isfield(f, 'tranMultR')
    trMult = f.tranMultR(p);
end
% Pressure dependent pore volume multiplier
if isfield(f, 'pvMultR')
    pvMult =  f.pvMultR(p);
    pvMult0 = f.pvMultR(p0);
end
if isfield(f, 'transMult')
   transMult = f.transMult(p); 
end
% Check for capillary pressure (p_cow)
pcOW = 0;
if isfield(f, 'pcOW')
    pcOW  = f.pcOW(sW);
end
% Gravity contribution, assert that it is aligned with z-dir
grav = gravity();
gdz = s.Grad(G.cells.centroids) * grav';

% Compute transmissibility
T = s.T.*transMult;

% Evaluate relative permeability
sO  = 1 - sW;
sO0 = 1 - sW0;

[krW, krO] = model.evaluteRelPerm({sW, sO});

% Multipliers due to polymer
mixpar = f.mixPar;
cbar   = c/f.cmax;
a = f.muWMult(f.cmax).^(1-mixpar);
b = 1./(1-cbar+cbar./a);
muWMult = b.*f.muWMult(c).^mixpar;
permRed = 1 + ((f.rrf-1)./f.adsMax).*effads(c, cmax, f);
muWMult  = muWMult.*permRed;

% Water props
bW     = f.bW(p);
rhoW   = bW.*f.rhoWS;
% rhoW on face, average of neighboring cells
rhoWf  = s.faceAvg(rhoW);
muW    = f.muW(p-pcOW);
muWeff = muWMult.*muW;
mobW   = trMult.*krW./muWeff;
dpW     = s.Grad(p-pcOW) - rhoWf.*gdz;
% water upstream-index
upcw = (double(dpW)<=0);
vW   = - s.faceUpstr(upcw, mobW).*T.*dpW;
bWvW = s.faceUpstr(upcw, bW).*vW;
if any(bW < 0)
    warning('Negative water compressibility present!')
end

% Polymer
mobP = (mobW.*c)./(a + (1-a)*cbar);
vP   = - s.faceUpstr(upcw, mobP).*s.T.*dpW;
bWvP = s.faceUpstr(upcw, bW).*vP;

% oil props
bO     = f.bO(p);
rhoO   = bO.*f.rhoOS;
rhoOf  = s.faceAvg(rhoO);
dpO    = s.Grad(p) - rhoOf.*gdz;
% oil upstream-index
upco = (double(dpO)<=0);
if isfield(f, 'BOxmuO')
    % mob0 is already multplied with b0
    mobO   = trMult.*krO./f.BOxmuO(p);
    bOvO   = - s.faceUpstr(upco, mobO).*T.*dpO;
    vO = bOvO./s.faceUpstr(upco, bO);
else
    mobO   = trMult.*krO./f.muO(p);
    vO = - s.faceUpstr(upco, mobO).*T.*dpO;
    bOvO   = s.faceUpstr(upco, bO).*vO;
end
if any(bO < 0)
    warning('Negative oil compressibility present!')
end

if model.outputFluxes
    state = model.storeFluxes(state, vW, vO, vP);
end

if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, []);
    state = model.storeMobilities(state, mobW, mobO, mobP);
    state = model.storeUpstreamIndices(state, upcw, upco, []);
end

% EQUATIONS ---------------------------------------------------------------
% oil:
eqs{1} = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*f.bO(p0).*sO0) + s.Div(bOvO);

% water:
eqs{2} = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*f.bW(p0).*sW0 ) + s.Div(bWvW);

% polymer in water:
poro =  s.pv./G.cells.volumes;
eqs{3} = (s.pv.*(1-f.dps)/dt).*(pvMult.*bW.*sW.*c - ...
   pvMult0.*f.bW(p0).*sW0.*c0) + (s.pv/dt).* ...
   ( f.rhoR.*((1-poro)./poro).*(effads(c, cmax, f)- ...
     effads(c0, cmax0, f)) ) + s.Div(bWvP);

names = {'oil', 'water', 'polymer'};
types = {'cell', 'cell', 'cell'};
% well equations
if ~isempty(W)
    if ~opt.reverseMode
        wc   = vertcat(W.cells);
        pw   = p(wc);
        rhos = [f.rhoWS, f.rhoOS];
        bw   = {bW(wc), bO(wc)};
        mw   = {mobW(wc), mobO(wc)};
        s    = {sW, 1 - sW};

        wm = WellModel();
        [cqs, weqs, ctrleqs, wc, state.wellSol] = ....
           wm.computeWellFlux(model, W, wellSol, pBH, {qWs, qOs}, ...
           pw, rhos, bw, mw, s, {}, 'nonlinearIteration', opt.iteration);
        eqs(4:5) = weqs;
        eqs{7} = ctrleqs;
        
        eqs{1}(wc) = eqs{1}(wc) - cqs{2};
        eqs{2}(wc) = eqs{2}(wc) - cqs{1};
        
        % Polymer well equations
        [~, wciPoly, iInxW] = getWellPolymer(W);
        cw        = c(wc);
        cw(iInxW) = wciPoly;
        cbarw     = cw/f.cmax;
        
        % Divide away water mobility and add in polymer
        bWqP = cw.*cqs{1}./(a + (1-a).*cbarw);
        eqs{3}(wc) = eqs{3}(wc) - bWqP;
        
        % Well polymer rate for each well is water rate in each perforation
        % multiplied with polymer concentration in that perforated cell.
        perf2well = getPerforationToWellMapping(W);
        eqs{6} = -qWPoly;
        for i = 1:numel(W)
            eqs{6}(i) = eqs{6}(i) + sum(cqs{1}(perf2well == i).*cw(perf2well == i));
        end
        
        names(4:7) = {'oilWells', 'waterWells', 'polymerWells', ...
           'closureWells'};
        types(4:7) = {'perf', 'perf', 'perf', 'well'};
    else
        % in reverse mode just gather zero-eqs of correct size
        for eqn = 4:7
            nw = numel(state0.wellSol);
            zw = double2ADI(zeros(nw,1), p0);
            eqs(4:7) = {zw, zw, zw, zw};
        end
        names(4:7) = {'empty', 'empty', 'empty', 'empty'};
        types(4:7) = {'none', 'none', 'none', 'none'};
    end
end
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
problem.iterationNo = opt.iteration;
end

%--------------------------------------------------------------------------


function [wPoly, wciPoly, iInxW] = getWellPolymer(W)
    if isempty(W)
        wPoly = [];
        wciPoly = [];
        iInxW = [];
        return
    end
    inj   = vertcat(W.sign)==1;
    polInj = cellfun(@(x)~isempty(x), {W(inj).poly});
    wPoly = zeros(nnz(inj), 1);
    wPoly(polInj) = vertcat(W(inj(polInj)).poly);
    wciPoly = rldecode(wPoly, cellfun(@numel, {W(inj).cells}));
    
    % Injection cells
    nPerf = cellfun(@numel, {W.cells})';
    nw    = numel(W);
    perf2well = rldecode((1:nw)', nPerf);
    compi = vertcat(W.compi);
    iInx  = rldecode(inj, nPerf);
    iInx  = find(iInx);
    iInxW = iInx(compi(perf2well(iInx),1)==1);
end

% Effective adsorption, depending of desorption or not
function y = effads(c, cmax, f)
   if f.adsInx == 2
      y = f.ads(max(c, cmax));
   else
      y = f.ads(c);
   end
end


