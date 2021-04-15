function [problem, state] = transportEquationCompositional(state0, state, model, dt, drivingForces, varargin)
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

opt = struct('Verbose',     mrstVerbose,...
            'reverseMode', false,...
            'computeFlash', true, ...
            'resOnly',     false,...
            'solveForWater', false, ...
            'iteration',   -1);

opt = merge_options(opt, varargin{:});

% Shorter names for some commonly used parts of the model and forces.
s = model.operators;
f = model.fluid;
W = drivingForces.W;
useMassFlux = state.massFlux;

fluid = model.fluid;
compFluid = model.EOSModel.fluid;

% if opt.computeFlash
    % Compute flash
%     state = model.computeFlash(state, dt, opt.iteration);
% end
state.eos.packed = model.EOSModel.getPropertiesFastAD(state.pressure, state.T, state.x, state.y);

% Properties at current timestep
[p, sW, z, wellSol] = model.getProps(state, ...
    'pressure', 'water', 'z', 'wellSol');

% L = state.L;
temp = state.T;

[p0, sW0, sO0, sG0, z0] = model.getProps(state0, ...
    'pressure', 'water', 'oil', 'gas', 'z');
z = expandMatrixToCell(z);
z0 = expandMatrixToCell(z0);


ncomp = numel(z);

if model.water
    [sW, z{1:ncomp-1}] = initVariablesADI(sW, z{1:ncomp-1});
    primaryVars = {'sW', compFluid.names{1:end-1}};
else
    [z{1:ncomp-1}] = initVariablesADI(z{1:ncomp-1});
    primaryVars = compFluid.names(1:end-1);
    sW = 0;
end

z{end} = 1;
for i = 1:(ncomp-1)
    z{end} = z{end} - z{i};
end
temp0 = state0.T;

[xM,  yM,  sO,  sG,  rhoO,  rhoG, muO, muG] = model.computeTwoPhaseFlowProps(state, p, temp, z);
[xM0, yM0, sO0, sG0, rhoO0, rhoG0] = model.computeTwoPhaseFlowProps(state0, p0, temp0, z0);


if model.water
    bW     = fluid.bW(p);
    bW0    = fluid.bW(p0);
    rhoW   = bW.*fluid.rhoWS;
    rhoW0 = bW0.*fluid.rhoWS;
    rhoT = sO.*rhoO + sG.*rhoG + sW.*rhoW;
    fsW = sW.*rhoW./rhoT;
else
    rhoT = sO.*rhoO + sG.*rhoG;
end
state.rhoT_transport = double(rhoT);

fsO = sO.*rhoO./rhoT;
fsG = sG.*rhoG./rhoT;

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

if model.water
    [krW, krO, krG] = model.evaluateRelPerm({sW, sO, sG});
    krW = mobMult.*krW;
else
    [krO, krG] = model.evaluateRelPerm({sO, sG});
end

krO = mobMult.*krO;
krG = mobMult.*krG;

% Compute transmissibility
T = s.T.*transMult;

% Gravity gradient per face
gdz = model.getGravityGradient();

% Water flux
if model.water
    muW = f.muW(p);

    rhoWf  = s.faceAvg(rhoW);
    mobW   = krW./muW;
    Gw = rhoWf.*gdz;
end

% muO = model.EOSModel.computeViscosity(p, rhoO, temp, x, true);
% muG = model.EOSModel.computeViscosity(p, rhoG, temp, y, false);

% Oil flux
rhoOf  = s.faceAvg(sO.*rhoO)./max(s.faceAvg(sO), 1e-8);
mobO   = krO./muO;

% Gas flux
rhoGf  = s.faceAvg(sG.*rhoG)./max(s.faceAvg(sG), 1e-8);
mobG   = krG./muG;

% splitting starts here
Go = rhoOf.*gdz;
Gg = rhoGf.*gdz;

flux = sum(state.flux, 2);
vT = flux(model.operators.internalConn);


% Stored upstream indices
if model.staticUpwind 
    flag = state.upstreamFlag;
else
    if useMassFlux
        if model.water
            gg = {Gw, Go, Gg};
            mg = {rhoW.*mobW, rhoO.*mobO, rhoG.*mobG};
        else
            gg = {Go, Gg};
            mg = {rhoO.*mobO, rhoG.*mobG};
        end
    else
        if model.water
            gg = {Gw, Go, Gg};
            mg = {mobW, mobO, mobG};
        else
            gg = {Go, Gg};
            mg = {mobO, mobG};
        end
    end
    flag = getSaturationUpwind(model.upwindType, state, gg, vT, s.T, mg, s.faceUpstr);
end
upco  = flag(:, model.water + 1);
upcg  = flag(:, model.water + 2);

if useMassFlux
    rhoT_p = state.rhoT;

    if model.water
        upcw  = flag(:, 1);
        rWmobWf = s.faceUpstr(upcw, mobW.*rhoW);
        rhoT0 = sO0.*rhoO0 + sG0.*rhoG0 + sW0.*rhoW0;
    else
        rWmobWf = 0;
        rhoT0 = sO0.*rhoO0 + sG0.*rhoG0;
    end
    rOmobOf = s.faceUpstr(upco, mobO.*rhoO);
    rGmobGf = s.faceUpstr(upcg, mobG.*rhoG);

    totRhoMob = rWmobWf + rOmobOf + rGmobGf;

    F_o = rOmobOf./totRhoMob;
    F_g = rGmobGf./totRhoMob;

    if model.water
        F_w = rWmobWf./totRhoMob;
        
        rOvO = F_o.*(vT + s.T.*rWmobWf.*(Go - Gw) + T.*rGmobGf.*(Go - Gg));
        rGvG = F_g.*(vT + s.T.*rWmobWf.*(Gg - Gw) + T.*rOmobOf.*(Gg - Go));
        rWvW = F_w.*(vT + s.T.*rOmobOf.*(Gw - Go) + T.*rGmobGf.*(Gw - Gg));
    else
        rOvO = F_o.*(vT + T.*rGmobGf.*(Go - Gg));
        rGvG = F_g.*(vT + T.*rOmobOf.*(Gg - Go));
    end
    fsO0 = sO0.*rhoO0./rhoT0;
    fsG0 = sG0.*rhoG0./rhoT0;
    
    if isfield(state0, 'rhoT')
        rhoT0_p = state0.rhoT;
    else
        rhoT0_p = state.rhoT0;
    end
else
    mobOf = s.faceUpstr(upco, mobO);
    mobGf = s.faceUpstr(upcg, mobG);
    
    rhoOf = s.faceUpstr(upco, rhoO);
    rhoGf = s.faceUpstr(upcg, rhoG);
    if model.water
        upcw  = flag(:, 1);
        mobWf = s.faceUpstr(upcw, mobW);
        rhoWf = s.faceUpstr(upcw, rhoW);
    else
        mobWf = 0;
        rhoWf = 0;
    end
    totMob = mobWf + mobOf + mobGf;
    F_o = mobOf./totMob;
    F_g = mobGf./totMob;

    if model.water
        F_w = mobWf./totMob;
        rOvO = rhoOf.*F_o.*(vT + s.T.*mobWf.*(Go - Gw) + T.*mobGf.*(Go - Gg));
        rGvG = rhoGf.*F_g.*(vT + s.T.*mobWf.*(Gg - Gw) + T.*mobOf.*(Gg - Go));
        rWvW = rhoWf.*F_w.*(vT + s.T.*mobOf.*(Gw - Go) + T.*mobGf.*(Gw - Gg));
    else
        rOvO = rhoOf.*F_o.*(vT + T.*mobGf.*(Go - Gg));
        rGvG = rhoGf.*F_g.*(vT + T.*mobOf.*(Gg - Go));
    end
end
[eqs, types, names] = deal(cell(1, ncomp-1+model.water));

eqOffset = 0;
if model.water
    
    if useMassFlux
        fsW0 = sW0.*rhoW0./rhoT0;
        eqs{1} = (s.pv/dt).*( rhoT_p.*pvMult.*fsW - rhoT0_p.*pvMult0.*fsW0 ) + s.Div(rWvW);
    else
        eqs{1} = (s.pv/dt).*( pvMult.*rhoW.*sW - pvMult0.*rhoW0.*sW0) + s.Div(rWvW);
    end
    names{1} = 'water';
    types{1} = 'cell';
    eqOffset = 1;
end

tmpeqs = cell(1, ncomp);
tmpnames = cell(1, ncomp);

for i = 1:ncomp%(ncomp - 1)
    if i < ncomp
        names{i + eqOffset} = compFluid.names{i};
        types{i + eqOffset} = 'cell';

        if useMassFlux
            eqs{i+eqOffset} = (s.pv/dt).*( ...
                            pvMult .*(rhoT_p .*(fsO .*xM {i} +  fsG.* yM{i} )) -...
                            pvMult0.*(rhoT0_p.*(fsO0.*xM0{i} + fsG0.* yM0{i}))...
                        )...
                  + s.Div(rOvO.*s.faceUpstr(upco, xM{i}) + rGvG.*s.faceUpstr(upcg, yM{i}));
        else
            eqs{i+eqOffset} = (s.pv/dt).*( ...
                        rhoO.*pvMult.*sO.*xM{i} - rhoO0.*pvMult0.*sO0.*xM0{i} + ...
                        rhoG.*pvMult.*sG.*yM{i} - rhoG0.*pvMult0.*sG0.*yM0{i} ...
                        )...
                  + s.Div(rOvO.*s.faceUpstr(upco, xM{i}) + rGvG.*s.faceUpstr(upcg, yM{i}));
        end
        if model.water
            pureWater = double(sW) == 1;
            if any(pureWater)
                % Cells with pure water should just retain their composition to
                % avoid singular systems
                eqs{i+eqOffset}(pureWater) = eqs{i+eqOffset}(pureWater) + ...
                                1e-3*(z{i}(pureWater) - double(z{i}(pureWater)));
            end
        end
    end
    
    tmpnames{i} = compFluid.names{i};

    tmpeqs{i+eqOffset} = (s.pv/dt).*( ...
                        rhoO.*pvMult.*sO.*xM{i} - rhoO0.*pvMult0.*sO0.*xM0{i} + ...
                        rhoG.*pvMult.*sG.*yM{i} - rhoG0.*pvMult0.*sG0.*yM0{i} ...
                        )...
                  + s.Div(rOvO.*s.faceUpstr(upco, xM{i}) + rGvG.*s.faceUpstr(upcg, yM{i}));

end

% Finally, add in and setup well equations
if ~isempty(W)
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
    
    if useMassFlux
        rOmobOw = rhoO(wc).*mobOw;
        rGmobGw = rhoG(wc).*mobGw;

        if model.water
            mobWw = mobW(wc);
            rWmobWw = rhoW(wc).*mobWw;
            totMobw = rWmobWw + rOmobOw + rGmobGw;
            f_w_w = rWmobWw./totMobw;
            f_w_w(isInj) = compPerf(isInj, 1);
            rWqW = f_w_w.*wflux;
            eqs{1}(wc) = eqs{1}(wc) - rWqW;
        else
            totMobw = rOmobOw + rGmobGw;
        end
        f_o_w = rOmobOw./totMobw;
        f_g_w = rGmobGw./totMobw;

        f_o_w(isInj) = compPerf(isInj, 1 + model.water);
        f_g_w(isInj) = compPerf(isInj, 2 + model.water);

        rOqO = f_o_w.*wflux;
        rGqG = f_g_w.*wflux;
    else
        mflux = sum(vertcat(wellSol.massFlux), 2);
        rhoOw = rhoO(wc);
        rhoGw = rhoG(wc);

        if model.water
            mobWw = mobW(wc);
            rhoWw = rhoW(wc);
            totMobw = mobWw + mobOw + mobGw;
            f_w_w = mobWw./totMobw;
            
            f_w_w(isInj) = compPerf(isInj, 1);
            rWqW = rhoWw.*f_w_w.*wflux;
            eqs{1}(wc) = eqs{1}(wc) - rWqW;
        else
            totMobw = mobOw + mobGw;
        end
        f_o_w = mobOw./totMobw;
        f_g_w = mobGw./totMobw;

        f_o_w(isInj) = compPerf(isInj, 1 + model.water);
        f_g_w(isInj) = compPerf(isInj, 2 + model.water);

        rOqO = rhoOw.*f_o_w.*wflux;
        rGqG = rhoGw.*f_g_w.*wflux;
        
        if 1
            rOqO(isInj) = compPerf(isInj, 1 + model.water).*mflux(isInj);
            rGqG(isInj) = compPerf(isInj, 2 + model.water).*mflux(isInj);
        end
    end

    sources = cell(ncomp, 1);
    compSrc = zeros(numel(wc), ncomp);
    for i = 1:ncomp
     
        ix = i + eqOffset;
        src =       (rOqO.*isInj + rGqG.*isInj).*w_comp(:, i)...
                   + rOqO.*(~isInj).*x_comp{i} + rGqG.*(~isInj).*y_comp{i};
        if i <= ncomp - 1
            eqs{ix}(wc) = eqs{ix}(wc) - src;
        end
        compSrc(:, i) = double(src);
        sources{i} = src;    
        tmpeqs{ix}(wc) = tmpeqs{ix}(wc) - src;
    end
    
    
    qg = double(rGqG)./fluid.rhoGS;
    qo = double(rOqO)./fluid.rhoOS;
    if model.water
        qw = double(rWqW)./fluid.rhoWS;
    end
    for i = 1:numel(W)
        state.wellSol(i).components = (compSrc(perf2well == i, :));
        state.wellSol(i).qGs = sum(qg(perf2well == i));
        state.wellSol(i).qOs = sum(qo(perf2well == i));
        if model.water
            state.wellSol(i).qWs = sum(qw(perf2well == i));
        end
    end
end

if model.water
    wscale = dt./(s.pv*mean(double(rhoW)));
    eqs{1} = eqs{1}.*wscale;
end

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

% state.s(:, 1 + model.water) = sL;
% state.s(:, 2 + model.water) = sV;

% Store pseudo-densities + saturations, accounting for volume error due to
% splitting of pressure and transport. Maybe move this into update
% functions?
if useMassFlux
    sL = double(sO);
    sV = double(sG);
    sA = double(sW);

    mG = double(rhoT_p).*double(fsG);
    mO = double(rhoT_p).*double(fsO);

    ro = mO./sL;
    bado = ~isfinite(ro);
    Ro = double(rhoO);
    ro(bado) = Ro(bado);

    rg = mG./sV;
    badg = ~isfinite(rg);
    Rg = double(rhoG);
    rg(badg) = Rg(badg);

    if model.water
        mW = double(rhoT_p).*double(fsW);

        rw = mW./sA;
        badw = ~isfinite(rw);
        Rw = double(rhoW);
        rw(badw) = Rw(badw);

        state.rho = [rw, ro, rg];
    else
        state.rho = [ro, rg];
    end
else
    if model.water
        state.rho = [double(rhoW), double(rhoO), double(rhoG)];
    else
        state.rho = [double(rhoO), double(rhoG)];
    end
end
problem = problem.assembleSystem();

if opt.iteration == 1
    state.residuals = [];
end
state.residuals = [state.residuals; cellfun(@(x) norm(double(x), inf), tmpeqs)];
state.residualNames = tmpnames;
end
