function [problem, state] = equationsOilWater(state0, state, model, dt, drivingForces, varargin)
% Get linearized problem for oil/water system with black oil-style
% properties
opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'iteration', -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.Wells;

% Operators, grid and fluid model.
s = model.operators;
G = model.G;
f = model.fluid;

% Properties at current timestep
[p, sW, wellSol] = model.getProps(state, 'pressure', 'water', 'wellsol');
% Properties at previous timestep
[p0, sW0] = model.getProps(state0, 'pressure', 'water');

pBH    = vertcat(wellSol.bhp);
qWs    = vertcat(wellSol.qWs);
qOs    = vertcat(wellSol.qOs);

% Initialize independent variables.
if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [p, sW, qWs, qOs, pBH] = ...
            initVariablesADI(p, sW, qWs, qOs, pBH);
    else
        [p0, sW0, tmp, tmp, tmp] = ...
            initVariablesADI(p0, sW0,          ...
            zeros(size(qWs)), ...
            zeros(size(qOs)), ...
            zeros(size(pBH)));                          %#ok
        clear tmp
    end
end
% We will solve for pressure, water saturation (oil saturation follows via
% the definition of saturations) and well rates + bhp.
primaryVars = {'pressure', 'sW', 'qWs', 'qOs', 'bhp'};

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
% Gravity contribution
gdz = model.getGravityGradient();

% Compute transmissibility
T = s.T.*transMult;

% Evaluate relative permeability
sO  = 1 - sW;
sO0 = 1 - sW0;

[krW, krO] = model.evaluteRelPerm({sW, sO});

% Water props
bW     = f.bW(p);
rhoW   = bW.*f.rhoWS;
% rhoW on face, average of neighboring cells
rhoWf  = s.faceAvg(rhoW);
mobW   = trMult.*krW./f.muW(p);
dpW     = s.Grad(p-pcOW) - rhoWf.*gdz;
% water upstream-index
upcw = (double(dpW)<=0);
vW   = - s.faceUpstr(upcw, mobW).*T.*dpW;
bWvW = s.faceUpstr(upcw, bW).*vW;
if any(bW < 0)
    warning('Negative water compressibility present!')
end

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
    state = model.storeFluxes(state, vW, vO, []);
end

if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, []);
    state = model.storeMobilities(state, mobW, mobO, []);
    state = model.storeUpstreamIndices(state, upcw, upco, []);
end

% EQUATIONS ---------------------------------------------------------------
% oil:
eqs{1} = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*f.bO(p0).*sO0 ) + s.Div(bOvO);

% water:
eqs{2} = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*f.bW(p0).*sW0 ) + s.Div(bWvW);

eqs([2, 1]) = addFluxesFromSourcesAndBC(model, eqs([2, 1]), ...
                                               {p - pcOW, p},...
                                               {rhoW,     rhoO},...
                                               {mobW,     mobO}, ...
                                               {bW, bO},  ...
                                               {sW, sO}, ...
                                               drivingForces);

names = {'oil', 'water'};
types = {'cell', 'cell'};
% well equations
if ~isempty(W)
    if ~opt.reverseMode
        wc    = vertcat(W.cells);
        pw   = p(wc);
        rhos = [f.rhoWS, f.rhoOS];
        bw   = {bW(wc), bO(wc)};
        mw   = {mobW(wc), mobO(wc)};
        s = {sW(wc), 1 - sW(wc)};

        wm = WellModel();
        [cqs, weqs, ctrleqs, wc, state.wellSol]  = wm.computeWellFlux(model, W, wellSol, ...
                                             pBH, {qWs, qOs}, pw, rhos, bw, mw, s, {},...
                                             'nonlinearIteration', opt.iteration);
        eqs(3:4) = weqs;
        eqs{5} = ctrleqs;
        
        eqs{1}(wc) = eqs{1}(wc) - cqs{2};
        eqs{2}(wc) = eqs{2}(wc) - cqs{1};
        
        names(3:5) = {'oilWells', 'waterWells', 'closureWells'};
        types(3:5) = {'perf', 'perf', 'well'};
    else
        % in reverse mode just gather zero-eqs of correct size
        for eqn = 3:5
            nw = numel(state0.wellSol);
            zw = double2ADI(zeros(nw,1), p0);
            eqs(3:5) = {zw, zw, zw};
        end
        names(3:5) = {'empty', 'empty', 'empty'};
        types(3:5) = {'none', 'none', 'none'};
    end
end



problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end
