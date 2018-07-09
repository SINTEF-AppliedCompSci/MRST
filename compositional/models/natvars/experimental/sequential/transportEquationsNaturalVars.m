function [problem, state] = transportEquationsNaturalVars(state0, state, model, dt, drivingForces, varargin)
opt = struct('Verbose',     mrstVerbose,...
            'reverseMode', false,...
            'resOnly',     false,...
            'iteration',   -1);

opt = merge_options(opt, varargin{:});


% Shorter names for some commonly used parts of the model and forces.
s = model.operators;
W = drivingForces.W;

fluid = model.fluid;
compFluid = model.EOSModel.fluid;

% Properties at current timestep
[p, sW, sO, sG, x, y, temp, wellSol] = model.getProps(state, ...
    'pressure', 'water', 'so', 'sg', 'x', 'y', 'T', 'wellSol');
assert(all(p>0), 'Pressure must be positive for compositional model');

[p0, sW0, sO0, sG0, x0, y0, temp0, wellSol0] = model.getProps(state0, ...
    'pressure', 'water', 'so', 'sg', 'x', 'y', 'T', 'wellSol');

[pureLiquid, pureVapor, twoPhase] = model.getFlag(state);
if 1
    stol = 1e-8;
    pureWater = sO + sG < stol;
    sO(~pureVapor & pureWater) = stol;
    sG(~pureLiquid & pureWater) = stol;

    [pureLiquid0, pureVapor0, twoPhase0] = model.getFlag(state0);
    pureWater0 = sO0 + sG0 < stol;
    sO0(~pureVapor0 & pureWater0) = stol;
    sG0(~pureLiquid0 & pureWater0) = stol;
end


% if isfield(state, 'timestep') && opt.iteration == 1
%     p = state.pressure_full;
%     dt_frac = dt/state.timestep;
%     state.pressure = p.*dt_frac + p0.*(1-dt_frac);
% end

z = state.components;
x(~twoPhase, :) = z(~twoPhase, :);
y(~twoPhase, :) = z(~twoPhase, :);
x = ensureMinimumFraction(x);
[y, z_tol] = ensureMinimumFraction(y);
x = expandMatrixToCell(x);
y = expandMatrixToCell(y);

ncomp = model.EOSModel.fluid.getNumberOfComponents();
[xnames, ynames, cnames] = deal(model.EOSModel.fluid.names);
for i = 1:ncomp
    xnames{i} = ['v_', cnames{i}];
    ynames{i} = ['w_', cnames{i}];
end
twoPhaseIx = find(twoPhase);

wtmp = ones(nnz(twoPhase), 1);
w = cell(ncomp, 1);
[w{:}] = deal(wtmp);

nc = model.G.cells.num;
for i = 1:(ncomp-1)
    w{i} = y{i}(twoPhase);
end
so = sO(twoPhase);

X = sG;
X(pureLiquid) = sO(pureLiquid);

if model.water
    if ~opt.resOnly
        [X, x{1:ncomp-1}, sW, w{1:ncomp-1}, so] = model.AutoDiffBackend.initVariablesAD(...
         X, x{1:ncomp-1}, sW, w{1:ncomp-1}, so);
    end
    primaryVars = {'sGsO', xnames{1:end-1}, 'satw', ynames{1:end-1}, 'sato'};

else
    if ~opt.resOnly
        [X, x{1:ncomp-1}, w{1:ncomp-1}, so] = model.AutoDiffBackend.initVariablesAD(...
         X, x{1:ncomp-1}, w{1:ncomp-1}, so);
    end
    primaryVars = {'sGsO', xnames{1:end-1}, ynames{1:end-1}, 'sato'};

end

sample = x{1};

sO = model.AutoDiffBackend.convertToAD(sO, sample);
sO(twoPhase) = so;
sO(pureLiquid) = X(pureLiquid);
sG = model.AutoDiffBackend.convertToAD(sG, sample);
sG(~pureLiquid) = X(~pureLiquid);

if model.water
    sT = sO + sG + sW;
    sT0 = sO0 + sG0 + sW0;
    sW = sW./sT;
    sW0 = sW0./sT0;
else
    sT = sO + sG;
    sT0 = sO0 + sG0;
end
sT = max(sT, 1e-8);
sT0 = max(sT0, 1e-8);

sO = sO./sT;
sG = sG./sT;

sO0 = sO0./sT0;
sG0 = sG0./sT0;
x{end} = 1;
w{end} = 1;
for i = 1:ncomp-1
    x{end} = x{end} - x{i};
    if any(twoPhase)
        w{end} = w{end}-w{i};
    end
end

for i = 1:ncomp
    y{i} = x{i};
    if any(twoPhase)
        y{i}(twoPhase) = w{i};
    end
    x{i}(pureVapor) = double(x{i}(pureVapor));
    y{i}(pureLiquid) = double(y{i}(pureLiquid));
end

cellJacMap = cell(numel(primaryVars), 1);
% toReorder = 1:nc;

if isempty(twoPhaseIx) || opt.resOnly
    reorder = [];
else
    n2ph = nnz(twoPhase);
    nVars = sum(sample.getNumVars());
    reorder = 1:nVars;
    start = twoPhaseIx + nc;
    stop = (nVars-n2ph+1):nVars;
    
    reorder(start) = stop;
    reorder(stop) = start;
    
    offset = ncomp+model.water;
    for i = 1:ncomp
        cellJacMap{i + offset} = twoPhaseIx;
    end
end


% Compute properties and fugacity
[xM,  yM,  rhoO,  rhoG,  muO,  muG, f_L, f_V, xM0, yM0, rhoO0, rhoG0] = ...
                  model.getTimestepPropertiesEoS(state, state0, p, temp, x, y, z, sO, sG, cellJacMap);

if model.water
    sat = {sW, sO, sG};
else
    sat = {sO, sG};
end

if model.water
    [krW, krO, krG] = model.evaluateRelPerm(sat);
else
    [krO, krG] = model.evaluateRelPerm(sat);
end

for i = 1:ncomp
    xM{i} = xM{i}.*sT;
    yM{i} = yM{i}.*sT;
    
    xM0{i} = xM0{i}.*sT0;
    yM0{i} = yM0{i}.*sT0;
end

% Compute transmissibility
T = s.T;

% Gravity gradient per face
gdz = model.getGravityGradient();
% rhoOf  = s.faceAvg(sat{1+model.water}.*rhoO)./max(s.faceAvg(sat{1+model.water}), 1e-8);
% rhoGf  = s.faceAvg(sat{2+model.water}.*rhoG)./max(s.faceAvg(sat{2+model.water}), 1e-8);

rhoOf  = s.faceAvg(sO.*rhoO.*sT)./max(s.faceAvg(sO.*sT), 1e-8);
rhoGf  = s.faceAvg(sG.*rhoG.*sT)./max(s.faceAvg(sG.*sT), 1e-8);


% Oil flux
mobO   = krO./muO;

% Gas flux
mobG   = krG./muG;


% splitting starts here
Go = rhoOf.*gdz;
Gg = rhoGf.*gdz;
if isfield(fluid, 'pcOG')
    Gg = Gg - s.Grad(fluid.pcOG(sG));
end

vT = sum(state.flux(model.operators.internalConn, :), 2);
if model.water
    pW = p;
    pW0 = p0;
    bW = fluid.bW(pW);
    rhoW = bW.*fluid.rhoWS;
    rhoW0 = fluid.bW(pW0).*fluid.rhoWS;
    
    rhoWf  = s.faceAvg(rhoW);
    muW = fluid.muW(pW);
    mobW   = krW./muW;
    Gw = rhoWf.*gdz;
    
    if isfield(fluid, 'pcOW')
        pcOW  = fluid.pcOW(sW);
        Gw = Gw + s.Grad(pcOW);
    end

    gg = {Gw, Go, Gg};
    mg = {mobW, mobO, mobG};
else
    [rhoW, rhoW0, mobW, bW] = deal([]);
    gg = {Go, Gg};
    mg = {mobO, mobG};
end
flag = getSaturationUpwind(model.upwindType, state, gg, vT, s.T, mg, s.faceUpstr);
upco  = flag(:, model.water + 1);
upcg  = flag(:, model.water + 2);


mobOf = s.faceUpstr(upco, mobO);
mobGf = s.faceUpstr(upcg, mobG);

rhoOf = s.faceUpstr(upco, rhoO);
rhoGf = s.faceUpstr(upcg, rhoG);

if model.water
    upcw  = flag(:, 1);
    mobWf = s.faceUpstr(upcw, mobW);
    rhoWf = s.faceUpstr(upcw, rhoW);
else
    upcw = [];
    mobWf = 0;
    rhoWf = 0;
end

if model.extraStateOutput
    bO = rhoO./fluid.rhoOS;
    bG = rhoG./fluid.rhoGS;

    state = model.storebfactors(state, bW, bO, bG);
    state = model.storeMobilities(state, mobW, mobO, mobG);
    state = model.storeUpstreamIndices(state, upcw, upco, upcg);
end
state = model.storeDensities(state, rhoW, rhoO, rhoG);

totMob = mobWf + mobOf + mobGf;
totMob = max(totMob, 1e-8);
F_o = mobOf./totMob;
F_g = mobGf./totMob;

if model.water
    F_w = mobWf./totMob;
    vO = F_o.*(vT + T.*mobWf.*(Go - Gw) + T.*mobGf.*(Go - Gg));
    vG = F_g.*(vT + T.*mobWf.*(Gg - Gw) + T.*mobOf.*(Gg - Go));
    vW = F_w.*(vT + T.*mobOf.*(Gw - Go) + T.*mobGf.*(Gw - Gg));

    rOvO = rhoOf.*vO;
    rGvG = rhoGf.*vG;
    rWvW = rhoWf.*vW;
else
    vO = F_o.*(vT + T.*mobGf.*(Go - Gg));
    vG = F_g.*(vT + T.*mobOf.*(Gg - Go));
    rOvO = rhoOf.*vO;
    rGvG = rhoGf.*vG;
    vW = [];
end

state = model.storeFluxes(state, vW, vO, vG);
% rOvO = s.faceUpstr(upco, sT).*rOvO;
% rGvG = s.faceUpstr(upcg, sT).*rGvG;

pv = model.operators.pv;
pv0 = pv;
if isfield(fluid, 'pvMultR')
    pv = pv.*fluid.pvMultR(p);
    pv0 = pv0.*fluid.pvMultR(p0);
end
compFlux = zeros(size(model.operators.N, 1), ncomp);

% water equation + n component equations
[eqs, types, names] = deal(cell(1, 2*ncomp + model.water));
for i = 1:ncomp
    names{i} = compFluid.names{i};
    types{i} = 'cell';
    
    vi = rOvO.*s.faceUpstr(upco, xM{i}) + rGvG.*s.faceUpstr(upcg, yM{i});
    eqs{i} = (1/dt).*( ...
                    pv.*rhoO.*sO.*xM{i} - pv0.*rhoO0.*sO0.*xM0{i} + ...
                    pv.*rhoG.*sG.*yM{i} - pv0.*rhoG0.*sG0.*yM0{i}) + s.Div(vi);
      
   compFlux(:, i) = double(vi);
end
state.componentFluxes = compFlux;
if model.water
    state.massFlux = [double(rWvW), double(rOvO), double(rGvG)];
else
    state.massFlux = [double(rOvO), double(rGvG)];
end

if model.water
    wix = ncomp+1;
    eqs{wix} = (1/dt).*(pv.*rhoW.*sW.*sT - pv0.*rhoW0.*sW0.*sT0) + s.Div(s.faceUpstr(upcw, sT).*rWvW);
    names{wix} = 'water';
    types{wix} = 'cell';
    
    rho = {rhoW, rhoO, rhoG};
    mob = {mobW, mobO, mobG};
    pressures = {p, p, p};
else
    rho = {rhoO, rhoG};
    mob = {mobO, mobG};
    pressures = {p, p};
end

comps = cellfun(@(x, y) {x, y}, xM, yM, 'UniformOutput', false);

[eqs, state] = model.addBoundaryConditionsAndSources(eqs, names, types, state, ...
                                                 pressures, sat, mob, rho, ...
                                                 {}, comps, ...
                                                 drivingForces);


% Finally, add in and setup well equations
if ~isempty(W)
    mflux = sum(vertcat(wellSol.components), 2);
    wflux = sum(vertcat(wellSol.flux), 2);

    perf2well = getPerforationToWellMapping(W);
    wc    = vertcat(W.cells);
    w_comp = vertcat(W.components);
    a = w_comp(perf2well, :).*repmat(compFluid.molarMass, numel(wc), 1);
    w_comp = bsxfun(@rdivide, a, sum(a, 2));

    isInj = wflux > 0;
    compWell = vertcat(W.compi);
    compPerf = compWell(perf2well, :);

    x_comp = cellfun(@(v) v(wc), xM, 'UniformOutput', false);
    y_comp = cellfun(@(v) v(wc), yM, 'UniformOutput', false);

    mobOw = mobO(wc);
    mobGw = mobG(wc);

    rhoOw = rhoO(wc);
    rhoGw = rhoG(wc);

    if model.water
        mobWw = mobW(wc);
        rhoWw = rhoW(wc);
        totMobw = mobWw + mobOw + mobGw;
        totMobw = max(totMobw, 1e-8);
        f_w_w = mobWw./totMobw;
        f_w_w = sT(wc).*f_w_w;

        f_w_w(isInj) = compPerf(isInj, 1);
        rWqW = rhoWw.*f_w_w.*wflux;
        eqs{wix}(wc) = eqs{wix}(wc) - rWqW;
    else
        totMobw = mobOw + mobGw;
        totMobw = max(totMobw, 1e-8);
    end
    f_o_w = mobOw./totMobw;
    f_g_w = mobGw./totMobw;
    rOqO = rhoOw.*f_o_w.*wflux;
    rGqG = rhoGw.*f_g_w.*wflux;


    rOqO(isInj) = compPerf(isInj, 1 + model.water).*mflux(isInj);
    rGqG(isInj) = compPerf(isInj, 2 + model.water).*mflux(isInj);

    sources = cell(ncomp, 1);
    compSrc = zeros(numel(wc), ncomp);
    for i = 1:ncomp
        src =       (rOqO + rGqG).*w_comp(:, i).*isInj...
                   +(rOqO.*x_comp{i} + rGqG.*y_comp{i}).*~isInj;
        eqs{i}(wc) = eqs{i}(wc) - src;
        compSrc(:, i) = double(src);
        sources{i} = src;    
    end

    qg = double(rGqG)./fluid.rhoGS;
    qo = double(rOqO)./fluid.rhoOS;

    isProd = ~isInj;
    qo(isProd) = qo(isProd).*sT(wc(isProd));
    qg(isProd) = qg(isProd).*sT(wc(isProd));
    if model.water
        qw = double(rWqW)./fluid.rhoWS;
        qw(isProd) = qw(isProd).*sT(wc(isProd));
    end
    for i = 1:numel(W)
        state.wellSol(i).components = compSrc(perf2well == i, :);
        state.wellSol(i).qGs = sum(qg(perf2well == i));
        state.wellSol(i).qOs = sum(qo(perf2well == i));
        if model.water
            state.wellSol(i).qWs = sum(qw(perf2well == i));
        end

        state.wellSol(i).qTr = sum(wflux(perf2well == i));
        state.wellSol(i).status = true;
        state.wellSol(i).type = W(i).type;
        state.wellSol(i).val = W(i).val;
    end
end

for i = 1:ncomp
    ix = i + ncomp + model.water;
    names{ix}= ['f_', compFluid.names{i}];
    types{ix} = 'fugacity';
    eqs{ix} = (f_L{i}(twoPhase) - f_V{i}(twoPhase))/barsa;
    absent = state.components(twoPhase, i) <= 10*z_tol;
    if model.water
        absent = absent | pureWater(twoPhase);
    end
    if any(absent) && isa(eqs{ix}, 'ADI')
        eqs{ix}.val(absent) = 0;
    end    
end
massT = model.getComponentScaling(state0);
scale = (dt./s.pv)./massT;
if model.water
    wscale = dt./(s.pv*mean(double(rhoW0)));
    eqs{wix} = eqs{wix}.*wscale;
end

for i = 1:ncomp
    eqs{i} = eqs{i}.*scale;
end


if model.reduceLinearSystem
    problem = ReducedLinearizedSystem(eqs, types, names, primaryVars, state, dt);
    problem.keepNum = model.G.cells.num*(ncomp+model.water);
    problem.reorder = reorder;
else
    problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end
end



%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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
