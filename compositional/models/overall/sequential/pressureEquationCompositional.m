function [problem, state] = pressureEquationCompositional(state0, state, model, dt, drivingForces, varargin)
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
            'resOnly',     false,...
            'stabilizeSaturation', true,...
            'staticWells',  false, ...
            'propsPressure', [], ...
            'computeFlash', true, ...
            'iteration',   -1);

opt = merge_options(opt, varargin{:});

% Shorter names for some commonly used parts of the model and forces.
s = model.operators;
f = model.fluid;
W = drivingForces.W;
fluid = model.fluid;

% Property pressure different from flow potential
p_prop = opt.propsPressure;
otherPropPressure = ~isempty(p_prop);

if ~otherPropPressure && opt.computeFlash
    state = model.computeFlash(state, dt, opt.iteration);
end
% Properties at current timestep
[p, sW, z, wellSol] = model.getProps(state, ...
    'pressure', 'water', 'z', 'wellSol');

[p0, sW0, sO0, sG0] = model.getProps(state0, ...
    'pressure', 'water', 'oil', 'gas');
z = expandMatrixToCell(z);

bhp    = vertcat(wellSol.bhp);
qWs    = vertcat(wellSol.qWs);
qOs    = vertcat(wellSol.qOs);
qGs    = vertcat(wellSol.qGs);
ncomp = numel(z);

if model.water
    [p, qWs, qOs, qGs, bhp] = initVariablesADI(p, qWs, qOs, qGs, bhp);
    primaryVars = {'pressure', 'qWs', 'qOs', 'qGs', 'bhp'};
else
    [p, qOs, qGs, bhp] = initVariablesADI(p, qOs, qGs, bhp);
    primaryVars = {'pressure', 'qOs', 'qGs', 'bhp'};
    sW = 0;
end
z{end} = 1;
for i = 1:(ncomp-1)
    z{end} = z{end} - z{i};
end
temp = state.T;

if otherPropPressure
    [x, y, L] = model.getProps(state, 'x', 'y', 'L');
    Z_L = state.Z_L;
    Z_V = state.Z_V;
else
    p_prop = p;
    [x, y, L] = model.EOSModel.getPhaseFractionAsADI(state, p_prop, temp, z);
    [Z_L, Z_V] = model.EOSModel.getCompressibility(state, p_prop, temp, x, y);
end




x0 = expandMatrixToCell(state0.x);
y0 = expandMatrixToCell(state0.y);
temp0 = state0.T;
rhoO = model.EOSModel.computeDensity(p_prop, x, Z_L, temp, true);
rhoG = model.EOSModel.computeDensity(p_prop, y, Z_V, temp, false);

if isfield(state0, 'rho') && 0
    % Precomputed values
    rhoO0 = state0.rho(:, 2);
    rhoG0 = state0.rho(:, 3);
else
    rhoO0 = model.EOSModel.computeDensity(p0, x0, state0.Z_L, temp0, true);
    rhoG0 = model.EOSModel.computeDensity(p0, y0, state0.Z_V, temp0, false);
end

% if isfield(state0, 'rho')
% %     rhoO = rhoO - rhoO0 + state0.rho(:, 1 + model.water);
% %     rhoG = rhoG - rhoG0 + state0.rho(:, 2 + model.water);
%     rhoO = (rhoO./rhoO0).*state0.rho(:, 1 + model.water);
%     rhoG = (rhoG./rhoG0).*state0.rho(:, 2 + model.water);
% end

if model.useMassFlux && opt.stabilizeSaturation && isfield(state0, 'LP')
    w = dt./(state0.splitting.dt + dt);
    if 1
        L = w*L + (1-w)*state0.LP;
    else
        wc = vertcat(W.cells);
        L(wc) = w*L(wc) + (1-w)*state0.LP(wc);
    end
end

sL = (1 - sW).*L.*Z_L./(L.*Z_L + (1-L).*Z_V);
sO = sL;

sG = 1 - sW - sO;


% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p_prop, p0);

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

if model.water
    % Water flux
    muW = f.muW(p_prop);
    bW     = fluid.bW(p_prop);
    rhoW   = bW.*fluid.rhoWS;
    rhoW0 = fluid.bW(p0).*fluid.rhoWS;

    rhoWf  = s.faceAvg(rhoW);
    mobW   = krW./muW;
    dpW    = s.Grad(p) - rhoWf.*gdz;
    upcw  = (double(dpW)<=0);
    vW = -s.faceUpstr(upcw, mobW).*T.*dpW;
    rWvW = s.faceUpstr(upcw, rhoW).*vW;
else
    rWvW = [];
    upcw = [];
    bW = zeros(model.G.cells.num, 1);
end
muO = model.EOSModel.computeViscosity(p_prop, rhoO, temp, x, true);
muG = model.EOSModel.computeViscosity(p_prop, rhoG, temp, y, false);

% Oil flux
rhoOf  = s.faceAvg(sO.*rhoO)./max(s.faceAvg(sO), 1e-8);
mobO   = krO./muO;
dpO    = s.Grad(p) - rhoOf.*gdz;
upco  = (double(dpO)<=0);
vO = -s.faceUpstr(upco, mobO).*T.*dpO;
% vO = -s.faceAvg(mobO).*T.*dpO;

% Gas flux
rhoGf  = s.faceAvg(sG.*rhoG)./max(s.faceAvg(sG), 1e-8);
mobG   = krG./muG;
dpG    = s.Grad(p) - rhoGf.*gdz;
upcg  = (double(dpG)<=0);
vG = -s.faceUpstr(upcg, mobG).*T.*dpG;
% vG = -s.faceAvg(mobG).*T.*dpG;

rOvO = s.faceUpstr(upco, rhoO).*vO;
rGvG = s.faceUpstr(upcg, rhoG).*vG;

state.massFlux = double(rOvO + rGvG);
state.volFlux = double(vO + vG);


if model.useMassFlux
    state = model.storeFluxes(state, rWvW, rOvO, rGvG);
    state.massFlux = true;
else
    state = model.storeFluxes(state, vW, vO, vG);
    state.massFlux = false;
end
state = model.storeUpstreamIndices(state, upcw, upco, upcg);
state.splitting.dt = dt;
state.splitting.sL_p = sL;
% EQUATIONS -----------------------------------------------------------

% Finally, add in and setup well equations
if ~isempty(W)
    wm = WellModel();
    
    if ~opt.reverseMode
        % Store cell wise well variables in cell arrays and send to ewll
        % model to get the fluxes and well control equations.
        wc    = vertcat(W.cells);
        if opt.staticWells
            q = vertcat(state.wellSol.flux);
            
            if model.water
                qW = q(:, 1)./fluid.rhoWS;
                qO = q(:, 2)./fluid.rhoOS;
                qG = q(:, 3)./fluid.rhoGS;

                cqs = {qW, qO, qG};
            else
                qO = q(:, 1)./fluid.rhoOS;
                qG = q(:, 2)./fluid.rhoGS;

                cqs = {qO, qG};
            end
            wnames = {};
            wtypes = {};
            weqs = {};
        else
            perf2well = getPerforationToWellMapping(W);
            pw    = p(wc);
            w_comp = vertcat(W.components);
            wx = cell(1, ncomp);
            for i = 1:ncomp
                wx{i} = w_comp(perf2well, i);
            end
            
            rhoOw = rhoO(wc);
            rhoGw = rhoG(wc);
            
            tmp = vertcat(W.sign);
            sgn = tmp(perf2well);
            
            isInj = sgn > 0;

            bOw = rhoOw./fluid.rhoOS;
            bGw = rhoGw./fluid.rhoGS;
            
%             double(rhoOw./rhoO_wellcell)
%             double(rhoGw./rhoG_wellcell)
            
            if model.water
                rhows = [f.rhoWS, f.rhoOS, f.rhoGS];
                bw    = {bW(wc), bOw, bGw};
                mw    = {mobW(wc), mobO(wc), mobG(wc)};
                sat   = {sW(wc), sO(wc), sG(wc)};
                rates = {qWs, qOs, qGs};
            else
                rhows = [f.rhoOS, f.rhoGS];
                bw    = {bOw, bGw};
                mw    = {mobO(wc), mobG(wc)};
                sat   = {sO(wc), sG(wc)};
                rates = {qOs, qGs};
            end

            [cqs, weqs, ctrleq, wc, state.wellSol, cqr]  = wm.computeWellFlux(model, W, wellSol, ...
                bhp, rates, pw, rhows, bw, mw, sat, {},...
                'nonlinearIteration', opt.iteration);
            weqs = {weqs{:}, ctrleq};
            if model.water
                wnames = {'waterWells', 'oilWells', 'gasWells', 'closureWells'};
                wtypes = {'perf', 'perf', 'perf', 'well'};
            else
                wnames = {'oilWells', 'gasWells', 'closureWells'};
                wtypes = {'perf', 'perf', 'well'};
            end
        end
        cqd = cellfun(@double, cqs, 'UniformOutput', false);
        if model.water
            qTot = fluid.rhoWS.*cqs{1} + fluid.rhoOS.*cqs{2} + fluid.rhoGS.*cqs{3};
            fluxt = [fluid.rhoWS.*cqd{1}, fluid.rhoOS.*cqd{2}, fluid.rhoGS.*cqd{3}];
        else
            qTot = fluid.rhoOS.*cqs{1} + fluid.rhoGS.*cqs{2};
            fluxt = [fluid.rhoOS.*cqd{1}, fluid.rhoGS.*cqd{2}];
        end
        
        if ~model.useMassFlux
            c = cellfun(@double, cqr, 'UniformOutput', false);
            fluxt = [c{:}];
        end
        
        for i = 1:numel(W)
            wp = wm.perf2well == i;
            state.wellSol(i).flux = fluxt(wp, :);
        end
    end
    
end

if model.water
    state.rho = [double(rhoW), double(rhoO), double(rhoG)];
    state.rhoP = state.rho;
    state.kr = [double(krW), double(krO), double(krG)];
    rhoT = sO.*rhoO + sG.*rhoG + sW.*rhoW;
%     rhoT = rhoT - rhoT_mod;
    rhoT0 = sO0.*rhoO0 + sG0.*rhoG0 + sW0.*rhoW0;
    vT = rWvW + rOvO + rGvG;
    
%     vT_prev = sum(state.flux(model.operators.internalConn, :), 2);
%     vT = (vW + vO + vG).*s.faceUpstr(upcg, rhoT);
else
    state.rho = [double(rhoO), double(rhoG)];
    state.rhoP = state.rho;
    state.kr = [double(krO), double(krG)];
    rhoT = sO.*rhoO + sG.*rhoG;
    rhoT0 = sO0.*rhoO0 + sG0.*rhoG0;
    vT = rOvO + rGvG;
end

if model.useMassFlux 
    if isfield(state0, 'rhoT')
        rhoT0 = state0.rhoT;
    end
    state.rhoT0 = double(rhoT0);
    state.rhoT = double(rhoT);

end

eq = (s.pv./dt).*(pvMult.*rhoT - pvMult0.*rhoT0) + s.Div(vT);
if ~isempty(W)
    eq(wc) = eq(wc) - qTot;
else
    weqs = {};
    wnames = {};
    wtypes = {};
end


if ~model.useMassFlux
    compFluid = model.EOSModel.fluid;
    [xM,  yM,  sO,  sG] = model.computeTwoPhaseFlowProps(state, p, temp, z);
    [xM0, yM0, sO0, sG0] = model.computeTwoPhaseFlowProps(state0, p0, temp0, state0.components);
    L_ix = model.water + 1;
    V_ix = L_ix + 1;
    
    w_comp = vertcat(W.components);
    perf2well = getPerforationToWellMapping(W);
    a = w_comp(perf2well, :).*repmat(compFluid.molarMass, numel(wc), 1);
    w_comp = bsxfun(@rdivide, a, sum(a, 2));
    
    x_comp = cellfun(@(v) v(wc), xM, 'UniformOutput', false);
    y_comp = cellfun(@(v) v(wc), yM, 'UniformOutput', false);
    injO = double(cqs{L_ix}) > 0;
    injG = double(cqs{V_ix}) > 0;
    D = getScaling(model, state);
    
    eq = 0;
    for i = 1:ncomp
        
        e = (s.pv/dt).*( ...
            rhoO.*pvMult.*sO.*xM{i} - rhoO0.*pvMult0.*sO0.*xM0{i} + ...
            rhoG.*pvMult.*sG.*yM{i} - rhoG0.*pvMult0.*sG0.*yM0{i}) ...
            + s.Div(rOvO.*s.faceUpstr(upco, xM{i}) + rGvG.*s.faceUpstr(upcg, yM{i}));
        
        src =       (fluid.rhoOS.*cqs{L_ix}.*injO + fluid.rhoGS.*cqs{V_ix}.*injG).*w_comp(wm.perf2well, i)...
            + fluid.rhoOS.*~injO.*x_comp{i}.*cqs{L_ix} + fluid.rhoGS.*~injG.*y_comp{i}.*cqs{V_ix};
        
        e(wc) = e(wc) - src;
        eq = eq + D(:, i).*e;
    end
end
scale = (dt./s.pv);
eq = eq.*scale;

eqs = {eq, weqs{:}}; %#ok
names = {'pressure', wnames{:}};
types = {'cell', wtypes{:}};

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end

function D = getScaling(model, state)
    [p, z, temp] = model.getProps(state, 'pressure', 'z', 'T');
    ncomp = numel(z);
    [p, z{1:ncomp-1}] = initVariablesADI(p, z{1:ncomp-1});
    z{end} = 1;
    for i = 1:ncomp-1
        z{end} = z{end} - z{i};
    end
    
    
    [xM,  yM,  sO,  sG,  rhoO,  rhoG, muO, muG] = model.computeTwoPhaseFlowProps(state, p, temp, z);
    C = cell(ncomp, 1);
    for i = 1:ncomp
        C{i} = sO.*rhoO.*xM{i} + sG.*rhoG.*yM{i};
    end
    D = zeros(model.G.cells.num, ncomp);
    for i = 1:model.G.cells.num
        J = zeros(ncomp, ncomp);
        for j = 1:ncomp
            for k = 1:ncomp
                J(j, k) = C{j}.jac{k}(i, i);
            end
        end
        b = zeros(ncomp, 1);
        b(1) = 1/barsa;
        w = (J')\b;
        D(i, :) = w;
    end
    % 
    % peq = 0; for i = 1:ncomp; peq = peq+C{i}.*D(:, i); end;
end

% function D = getScaling(model, state)
%     [p, z, temp] = model.getProps(state, 'pressure', 'z', 'T');
%     sumVals = false;
%     if sumVals;
%         [z{1:end-1}] = initVariablesADI(z{1:end-1});
%         z{end} = 1;
%         for i = 1:numel(z)-1
%             z{end} = z{end} - z{i};
%         end
%     else
%         [z{:}] = initVariablesADI(z{:});
%     end
%     [xM,  yM,  sO,  sG,  rhoO,  rhoG, muO, muG] = model.computeTwoPhaseFlowProps(state, p, temp, z);
%     D = zeros(model.G.cells.num, numel(z));
%     
%     
%     for i = 1:numel(z)-sumVals
%         R = rhoO.*sO.*xM{i} + rhoG.*sG.*yM{i};
% 
% %         presentLiquid = double(sG) > 0;
% %         presentVapor = double(sO) > 0;
%         
% %         R = rhoO.*sO.*xM{i}.*presentLiquid + rhoG.*sG.*yM{i}.*presentVapor;
%         D(:, i) = 1./diag(R.jac{i});
%     end
%     
%     if sumVals
%         D(:, end) = 1 - sum(D(:, 1:end-1), 2);
%     end
% end










% mod = ThreePhaseCompositionalModel(model.G, model.rock, model.fluid, model.EOSModel.fluid);
% 
% probFIMP = mod.getEquations(state0, state, dt, drivingForces, varargin{:});
% 
% fimpEq = 0;
% 
% for i = 1:ncomp
%     fimpEq = fimpEq + probFIMP.equations{i}.*D(:, i);
% end
% fimpEq = fimpEq.*scale;
% 
% fimpEq.val./eq.val
