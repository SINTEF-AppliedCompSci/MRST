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
mixture = model.EOSModel.CompositionalMixture;

% Properties at current timestep
[sT, p, sW, z, temp, wellSol] = model.getProps(state, ...
    'sT', 'pressure', 'water', 'z', 'T', 'wellSol');
[p0, sW0, sO0, sG0, z0, temp0] = model.getProps(state0, ...
    'pressure', 'water', 'oil', 'gas', 'z', 'T');

z = ensureMinimumFraction(z);
z0 = ensureMinimumFraction(z0);

z = expandMatrixToCell(z);
z0 = expandMatrixToCell(z0);

ncomp = numel(z);
cnames = model.EOSModel.getComponentNames();

init = @(varargin) model.AutoDiffBackend.initVariablesAD(varargin{:});
if model.water
    [sT, z{1:ncomp-1}, sW] = init(sT, z{1:ncomp-1}, sW);
    primaryVars = {'sT', cnames{1:end-1}, 'sW'};
else
    [sT, z{1:ncomp-1}] = init(sT, z{1:ncomp-1});
    primaryVars = {'sT', cnames{1:end-1}};
end
z{end} = 1;
for i = 1:(ncomp-1)
    z{end} = z{end} - z{i};
end


[xM,  yM,  sO,  sG,  rhoO,  rhoG, muO, muG] = model.computeTwoPhaseFlowProps(state, p, temp, z);
[xM0, yM0, ~, ~, rhoO0, rhoG0]          = model.computeTwoPhaseFlowProps(state0, p0, temp0, z0);

[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

if model.water
    sO = (1-sW).*sO;
    sG = (1-sW).*sG;
    sat = {sW, sO, sG};
else
    sat = {sO, sG};
end

if model.water
    [krW, krO, krG] = model.evaluateRelPerm(sat);
    krW = mobMult.*krW;
else
    [krO, krG] = model.evaluateRelPerm(sat);
end

krO = mobMult.*krO;
krG = mobMult.*krG;

for i = 1:ncomp
    xM{i} = xM{i}.*sT;
    yM{i} = yM{i}.*sT;
end
% Compute transmissibility
T = transMult.*s.T;

% Gravity gradient per face
gdz = model.getGravityGradient();
rhoOf  = s.faceAvg(sO.*rhoO)./max(s.faceAvg(sO), 1e-8);
rhoGf  = s.faceAvg(sG.*rhoG)./max(s.faceAvg(sG), 1e-8);
% Phase mobility
mobO   = krO./muO;
mobG   = krG./muG;

Go = rhoOf.*gdz;
Gg = rhoGf.*gdz;

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
    sat{1} = sat{1}.*sT;
    gg = {Gw, Go, Gg};
    mob = {mobW, mobO, mobG};
    rho = {rhoW, rhoO, rhoG};
    pressures = {pW, p, p};

else
    [rhoW, rhoW0, mobW, bW] = deal([]);
    gg = {Go, Gg};
    mob = {mobO, mobG};
    rho = {rhoO, rhoG};
    pressures = {p, p};
end

if model.extraStateOutput
    bO = rhoO./fluid.rhoOS;
    bG = rhoG./fluid.rhoGS;

    state = model.storebfactors(state, bW, bO, bG);
    state = model.storeMobilities(state, mobW, mobO, mobG);
end
state = model.storeDensities(state, rhoW, rhoO, rhoG);

components = getComponentsTwoPhaseSimpleWater(model, rho, sT, xM, yM);

upstr = model.operators.faceUpstr;
[q_phase, q_components] = computeSequentialFluxes(...
    state, gg, vT, T, mob, rho, components, upstr, model.upwindType);


pv = pvMult.*model.operators.pv;
pv0 = pvMult0.*model.operators.pv;

% water equation + n component equations
[eqs, types, names] = deal(cell(1, ncomp + model.water));
for i = 1:ncomp
    names{i} = mixture.names{i};
    types{i} = 'cell';
    eqs{i} = (1/dt).*(pv.*rhoO.*sO.*xM{i} - pv0.*rhoO0.*sO0.*xM0{i} + ...
                      pv.*rhoG.*sG.*yM{i} - pv0.*rhoG0.*sG0.*yM0{i});
    if model.water
        pureWater = value(sG) == 0 & value(sO) == 0;
        if any(pureWater)
            % Cells with pure water should just retain their composition to
            % avoid singular systems
            eqs{i}(pureWater) = eqs{i}(pureWater) + ...
                            1e-3*(xM{i}(pureWater) - value(xM{i}(pureWater)));
        end
    end
end
if model.water
    wix = ncomp+1;
    eqs{wix} = (1/dt).*(pv.*rhoW.*sW.*sT - pv0.*rhoW0.*sW0);
    names{wix} = 'water';
    types{wix} = 'cell';
    state = model.storeFluxes(state, q_phase{:});
else
    state = model.storeFluxes(state, [], q_phase{:});
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
    a = w_comp(perf2well, :).*repmat(mixture.molarMass, numel(wc), 1);
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
        compSrc(:, i) = value(src);
        sources{i} = src;    
    end

    qg = value(rGqG)./fluid.rhoGS;
    qo = value(rOqO)./fluid.rhoOS;

    isProd = ~isInj;
    sTw = value(sT(wc(isProd)));
    qo(isProd) = qo(isProd).*sTw;
    qg(isProd) = qg(isProd).*sTw;
    if model.water
        qw = value(rWqW)./fluid.rhoWS;
        qw(isProd) = qw(isProd).*sTw;
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
% Apply scaling and assemble transport equations
for i = 1:(ncomp + model.water)
    vi = q_components{i};
    eqs{i} = s.AccDiv(eqs{i}, vi);
end

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end



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
