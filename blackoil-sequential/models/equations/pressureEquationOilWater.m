function [problem, state] = pressureEquationOilWater(state0, state, model, dt, drivingForces, varargin)

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'scaling', [],...
             'resOnly', false,...
             'history', [],...
             'iteration', -1, ...
             'stepOptions', []);  % Compatibility only

opt = merge_options(opt, varargin{:});

W = drivingForces.Wells;
perf2well = getPerforationToWellMapping(W);

assert(isempty(drivingForces.bc) && isempty(drivingForces.src))

s = model.operators;
f = model.fluid;
G = model.G;

[p, sW, wellSol] = model.getProps(state, 'pressure', 'water', 'wellsol');

[p0, sW0] = model.getProps(state0, 'pressure', 'water');


pBH    = vertcat(wellSol.bhp);
qWs    = vertcat(wellSol.qWs);
qOs    = vertcat(wellSol.qOs);

%Initialization of independent variables ----------------------------------

if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [p, qWs, qOs, pBH] = ...
            initVariablesADI(p, qWs, qOs, pBH);
    else
        assert(0, 'Backwards solver not supported for splitting');
    end
end
primaryVars = {'pressure', 'qWs', 'qOs', 'bhp'};

clear tmp
grav  = gravity;


%check for p-dependent porv mult:
pvMult = 1; pvMult0 = 1;
if isfield(f, 'pvMultR')
    pvMult =  f.pvMultR(p);
    pvMult0 = f.pvMultR(p0);
end

if 0 && isfield(wellSol, 'flux')
    % Linearize saturations in well cells to get mobilities at end of time
    % integration sort-of-right.
    bW = f.bW(p);
    bO = f.bO(p);
    
    flux = vertcat(wellSol.flux);
    wc = vertcat(W.cells);
    
    perfcells = wc(perf2well);    
    in = flux*dt./repmat(s.pv(perfcells), 1, 2);
    
    sW(perfcells) = sW(perfcells) + in(:, 1)./bW(perfcells) - in(:, 2)./bO(perfcells);
    sW = min(sW, 1);
    sW = max(sW, 0);
end
% -------------------------------------------------------------------------
sO  = 1 - sW;
sO0 = 1 - sW0;

[krW, krO] = model.evaluteRelPerm({sW, sO});

%dZ = s.grad(G.cells.centroids(:,3));
gdz = s.Grad(G.cells.centroids) * grav';

% Water
[bW, rhoW, mobW, dpW] = propsOW_water(sW, krW, gdz, f, p, s);

% water upstream-index
upcw = (double(dpW)<=0);
vW = - s.faceUpstr(upcw, mobW).*s.T.*dpW;
bWvW = s.faceUpstr(upcw, bW).*vW;


[bO, rhoO, mobO, dpO] = propsOW_oil(1 - sW, krO, gdz, f, p, s);
% oil upstream-index
upco = (double(dpO)<=0);
vO = - s.faceUpstr(upco, mobO).*s.T.*dpO;
bOvO = s.faceUpstr(upco, bO).*vO;

% These are needed in transport solver, so we output them regardless of
% any flags set in the model.
state = model.storeFluxes(state, vW, vO, []);
state = model.storeUpstreamIndices(state, upcw, upco, []);

% EQUATIONS ---------------------------------------------------------------
% oil:
oil = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*f.bO(p0).*sO0) + s.Div(bOvO);

% water:
wat = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*f.bW(p0).*sW0 ) + s.Div(bWvW);

[eqs, names, types] = deal({});

% well equations
if ~isempty(W)
    wc    = vertcat(W.cells);
    pw   = p(wc);
    rhos = [f.rhoWS, f.rhoOS];
    bw   = {bW(wc), bO(wc)};
    mw   = {mobW(wc), mobO(wc)};
    s = {sW, 1 - sW};

    wm = WellModel();
    [cqs, weqs, ctrleqs, wc, state.wellSol, cqr]  = wm.computeWellFlux(model, W, wellSol, ...
                                         pBH, {qWs, qOs}, pw, rhos, bw, mw, s, {},...
                                         'nonlinearIteration', opt.iteration);
    eqs(2:3) = weqs;
    eqs{4} = ctrleqs;
    
    qW = cqr{1};
    qO = cqr{2};
    
    oil(wc) = oil(wc) - cqs{2};
    wat(wc) = wat(wc) - cqs{1};

    names(2:4) = {'oilWells', 'waterWells', 'closureWells'};
    types(2:4) = {'perf', 'perf', 'well'};
end

eqs{1} = oil./bO + wat./bW;
names{1} = 'pressure';
types{1} = 'cell';

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

for i = 1:numel(W)
    wp = perf2well == i;
    state.wellSol(i).flux = [double(qW(wp)), double(qO(wp))];
end
end
