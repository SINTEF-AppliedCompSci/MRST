function [problem, state] = transportEquationOverallComposition(state0, state, model, dt, drivingForces, varargin)
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
[sT, p, sW, z, temp, wellSol] = model.getProps(state, ...
    'sT', 'pressure', 'water', 'z', 'T', 'wellSol');

[sT0, p0, sW0, sO0, sG0, z0, temp0, wellSol0] = model.getProps(state0, ...
    'sT', 'pressure', 'water', 'oil', 'gas', 'z', 'T', 'wellSol');
z = expandMatrixToCell(z);
z0 = expandMatrixToCell(z0);

ncomp = numel(z);
cnames = model.EOSModel.fluid.names;

if isfield(state, 'timestep') && opt.iteration == 1
    p = state.pressure_full;
    dt_frac = dt/state.timestep;
    state.pressure = p.*dt_frac + p0.*(1-dt_frac);
end

if model.EOSModel.fastDerivatives
    state.eos.packed = model.EOSModel.getPropertiesFastAD(state.pressure, state.T, state.x, state.y, state.components);
else
    state.eos.packed = struct();
end


if model.water
    [sT, z{1:ncomp-1}, sW] = initVariablesADI(sT, z{1:ncomp-1}, sW);
    primaryVars = {'sT', cnames{1:end-1}, 'sW'};
else
    [z{1:ncomp-1}, sT] = initVariablesADI(z{1:ncomp-1}, sT);
    primaryVars = {cnames{1:end-1}, 'sT'};
end


[xM,  yM,  sO,  sG,  rhoO,  rhoG, muO, muG] = model.computeTwoPhaseFlowProps(state, p, temp, z);
[xM0, yM0, ~, ~, rhoO0, rhoG0]          = model.computeTwoPhaseFlowProps(state0, p0, temp0, z0);

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
    
%     xM0{i} = xM0{i}.*sT0;
%     yM0{i} = yM0{i}.*sT0;
end


% disp([min(double(sT)), max(double(sT))])
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

vT = sum(state.flux(model.operators.internalConn, :), 2);
if model.water
    bW = model.fluid.bW(p);
    rhoW = bW.*model.fluid.rhoWS;
    rhoW0 = model.fluid.bW(p0).*model.fluid.rhoWS;
    
    rhoWf  = s.faceAvg(rhoW);
    muW = model.fluid.muW(p);
    mobW   = krW./muW;
    Gw = rhoWf.*gdz;

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
end

% rOvO = s.faceUpstr(upco, sT).*rOvO;
% rGvG = s.faceUpstr(upcg, sT).*rGvG;

pv = model.operators.pv;
pv0 = pv;
if isfield(fluid, 'pvMultR')
    pv = pv.*fluid.pvMultR(p);
    pv0 = pv0.*fluid.pvMultR(p0);
end


% water equation + n component equations
[eqs, types, names] = deal(cell(1, ncomp + model.water));
for i = 1:ncomp
    names{i} = compFluid.names{i};
    types{i} = 'cell';
      
    eqs{i} = (1/dt).*( ...
                    pv.*rhoO.*sO.*xM{i} - pv0.*rhoO0.*sO0.*xM0{i} + ...
                    pv.*rhoG.*sG.*yM{i} - pv0.*rhoG0.*sG0.*yM0{i}) ...
          + s.Div(rOvO.*s.faceUpstr(upco, xM{i}) + rGvG.*s.faceUpstr(upcg, yM{i}));
      
    if model.water
        pureWater = double(sG) == 0 & double(sO) == 0;
        if any(pureWater)
            % Cells with pure water should just retain their composition to
            % avoid singular systems
            eqs{i}(pureWater) = eqs{i}(pureWater) + ...
                            1e-3*(xM{i}(pureWater) - double(xM{i}(pureWater)));
        end
    end
end
if model.water
    wix = ncomp+1;
    eqs{wix} = (1/dt).*(pv.*rhoW.*sW.*sT - pv0.*rhoW0.*sW0) + s.Div(s.faceUpstr(upcw, sT).*rWvW);
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

if model.water
    wscale = dt./(s.pv*mean(double(rhoW0)));
    eqs{wix} = eqs{wix}.*wscale;
    
    scale = (dt./s.pv)./mean(double(sO0).*double(rhoO0) + double(sG0).*double(rhoG0) + double(sW0).*double(rhoW0));
else
    scale = (dt./s.pv)./mean(double(sO0).*double(rhoO0) + double(sG0).*double(rhoG0));
end
for i = 1:ncomp
    eqs{i} = eqs{i}.*scale;
end


problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
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
