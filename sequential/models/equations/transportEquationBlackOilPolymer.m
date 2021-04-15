function [problem, state] = transportEquationBlackOilPolymer(state0, state, model, dt, drivingForces, varargin)

opt = struct('Verbose'      , mrstVerbose, ...
             'reverseMode'  , false      , ...
             'scaling'      , []         , ...
             'resOnly'      , false      , ...
             'history'      , []         , ...
             'solveForWater', false      , ...
             'solveForOil'  , true       , ...
             'solveForGas'  , true       , ...
             'iteration'    , -1         , ...
             'stepOptions'  , []           );  % Compatibility only

opt = merge_options(opt, varargin{:});
assert(isempty(drivingForces.bc) && isempty(drivingForces.src))

W = drivingForces.W;
op = model.operators;
fluid = model.fluid;

disgas = model.disgas;
vapoil = model.vapoil;

% Polymer shear thinning/thickening
usingShear = isfield(fluid, 'plyshearMult');

assert(~(opt.solveForWater && opt.solveForOil && opt.solveForGas));

% Properties at current timestep
[p, sW, sG, c, cmax, rs, rv, wellSol] = model.getProps(state, ...
    'pressure', 'water', 'gas', 'polymer', 'polymermax', 'rs', 'rv', 'wellSol');

% Properties at previous timestep
[p0, sW0, sG0, c0, cmax0, rs0, rv0] = model.getProps(state0, ...
    'pressure', 'water', 'gas', 'polymer', 'polymermax', 'rs', 'rv');

% If timestep has been split relative to pressure, linearly interpolate in
% pressure.
pFlow = p;
if isfield(state, 'timestep')
    dt_frac = dt/state.timestep;
    p = p.*dt_frac + p0.*(1-dt_frac);
end

%Initialization of primary variables ----------------------------------
st  = model.getCellStatusVO(state,  1-sW-sG,   sW,  sG);
st0 = model.getCellStatusVO(state0, 1-sW0-sG0, sW0, sG0);
if ~opt.resOnly,
    if ~opt.reverseMode,
        % define primary varible x and initialize
        x = st{1}.*rs + st{2}.*rv + st{3}.*sG;

        [sW, x, c] = model.AutoDiffBackend.initVariablesAD(sW, x, c);

        % define sG, rs and rv in terms of x
        sG = st{2}.*(1-sW) + st{3}.*x;

        if disgas
            rsSat = fluid.rsSat(p);
            rs = (~st{1}).*rsSat + st{1}.*x;
        end
        if vapoil
            rvSat = fluid.rvSat(p);
            rv = (~st{2}).*rvSat + st{2}.*x;
        end
    else
        assert(0, 'Backwards solver not supported for splitting');
    end
end
if disgas || vapoil
    gvar = 'x';
else
    gvar = 'sG';
end
primaryVars = {'sW', gvar, 'polymer'};

% Evaluate relative permeability
sO  = 1 - sW  - sG;
sO0 = 1 - sW0 - sG0;
[krW, krO, krG] = model.evaluateRelPerm({sW, sO, sG});

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

% Modifiy relperm by mobility multiplier (if any)
krW = mobMult.*krW; krO = mobMult.*krO; krG = mobMult.*krG;

% Compute transmissibility
T = op.T.*transMult;

% Gravity gradient per face
gdz = model.getGravityGradient();

% Evaluate water and polymer properties
ads  = effads(c, cmax, model);
ads0 = effads(c0, cmax0, model);
[vW, vP, bW, muWMult, mobW, mobP, rhoW, pW, upcw, a, dpW] = ...
    getFluxAndPropsWaterPolymer_BO(model, p, sW, c, ads, krW, T, ...
    gdz, 'shear', false); % shear effect is not used in transport
bW0 = fluid.bW(p0);

% Evaluate oil properties
[vO, bO, mobO, rhoO, pO, upco, dpO] = getFluxAndPropsOil_BO(model, p, sO, krO, T, gdz, rs, ~st{1});
bO0 = getbO_BO(model, p0, rs0, ~st0{1});

% Evaluate gas properties
[vG, bG, mobG, rhoG, pG, upcg, dpG] = getFluxAndPropsGas_BO(model, p, sG, krG, T, gdz, rv, ~st{2});
bG0 = getbG_BO(model, p0, rv0, ~st0{2});

% Get total flux from state
flux = sum(state.flux, 2);
vT = flux(model.operators.internalConn);

% Sat dependent pressure terms
gp = op.Grad(p);
Gw = gp - dpW;
Go = gp - dpO;
Gg = gp - dpG;

[flag_v, flag_g] = getSaturationUpwind(model.upwindType, state, {Gw, Go, Gg}, vT, op.T, {mobW, mobO, mobG}, op.faceUpstr);

upcw  = flag_v(:, 1);
upco  = flag_v(:, 2);
upcg  = flag_v(:, 3);

% Upstream weighted face mobilities
mobWf = op.faceUpstr(upcw, mobW);
mobOf = op.faceUpstr(upco, mobO);
mobGf = op.faceUpstr(upcg, mobG);
mobPf = op.faceUpstr(upcw, c.*mobP);
if usingShear
    % The shear multipliers from the pressure solver are used
    mobWf = mobWf.*state.shearMult;
    mobPf = mobPf.*state.shearMult;
end
% Tot mob
totMob = mobOf + mobWf + mobGf;
totMob = max(totMob, sqrt(eps));

f_w = mobWf./totMob;
f_o = mobOf./totMob;
f_g = mobGf./totMob;
f_p = mobPf./totMob;

vW = f_w.*(vT + op.T.*mobOf.*(Gw - Go) + op.T.*mobGf.*(Gw - Gg));
vO = f_o.*(vT + op.T.*mobWf.*(Go - Gw) + op.T.*mobGf.*(Go - Gg));
vG = f_g.*(vT + op.T.*mobWf.*(Gg - Gw) + op.T.*mobOf.*(Gg - Go));
vP = f_p.*(vT + op.T.*mobOf.*(Gw - Go) + op.T.*mobGf.*(Gw - Gg));

bWvW = op.faceUpstr(upcw, bW).*vW;
bOvO = op.faceUpstr(upco, bO).*vO;
bGvG = op.faceUpstr(upcg, bG).*vG;
bWvP = op.faceUpstr(upcw, bW).*vP;

if disgas
    rsbOvO = op.faceUpstr(upco, rs).*bOvO;
end

if vapoil
    rvbGvG = op.faceUpstr(upcg, rv).*bGvG;
end

if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, bG);
    state = model.storeMobilities(state, mobW, mobO, mobG, mobP);
end
% well equations
if ~isempty(W)
    wflux = sum(vertcat(wellSol.flux), 2);
    cqs = vertcat(wellSol.cqs);

    perf2well = getPerforationToWellMapping(W);
    wc    = vertcat(W.cells);

    mobWw = mobW(wc);
    mobOw = mobO(wc);
    mobGw = mobG(wc);
    if usingShear
        % The shear multipliers from the pressure solver are used
        shearMultW = vertcat(state.wellSol.shearMult);
        mobWw = mobWw.*shearMultW;
    end
    totMobw = mobWw + mobOw + mobGw;


    isInj = wflux > 0;
    compWell = vertcat(W.compi);
    compPerf = compWell(perf2well, :);

    f_w_w = mobWw./totMobw;
    f_o_w = mobOw./totMobw;
    f_g_w = mobGw./totMobw;


    f_w_w(isInj) = compPerf(isInj, 1);
    f_o_w(isInj) = compPerf(isInj, 2);
    f_g_w(isInj) = compPerf(isInj, 3);


    bWqW = bW(wc).*f_w_w.*wflux;
    bOqO = bO(wc).*f_o_w.*wflux;
    bGqG = bG(wc).*f_g_w.*wflux;

    % Polymer well equations
    [~, wciPoly, iInxW] = getWellPolymer(W);
    cw        = c(wc);
    cw(iInxW) = wciPoly;
    bWqP      = cw.*bWqW;
    if 0
        bWqW(isInj) = cqs(isInj, 1);
        bOqO(isInj) = cqs(isInj, 2);
        bGqG(isInj) = cqs(isInj, 3);
        bWqP(isInj) = cqs(isInj, 4);
    end
    % Store well fluxes
    wflux_O = bOqO;
    wflux_W = bWqW;
    wflux_G = bGqG;
    wflux_P = bWqP;
    
    if disgas
        wflux_G = wflux_G + ~isInj.*bOqO.*rs(wc);
    end

    if vapoil
        wflux_O = wflux_O + ~isInj.*bGqG.*rv(wc);
    end

    for i = 1:numel(W)
        perfind = perf2well == i;
        state.wellSol(i).qOs = sum(value(wflux_O(perfind)));
        state.wellSol(i).qWs = sum(value(wflux_W(perfind)));
        state.wellSol(i).qGs = sum(value(wflux_G(perfind)));
        state.wellSol(i).qPs = sum(value(wflux_P(perfind)));
    end
end

[eqs, names, types] = deal(cell(1,3));
ix = 1;
if opt.solveForWater
    % water eq:
    water = (op.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + op.Div(bWvW);
    if ~isempty(W)
        water(wc) = water(wc) - wflux_W;
    end
    eqs{ix} = water;
    names{ix} = 'water';
    types{ix} = 'cell';
    ix = ix + 1;
end

if opt.solveForOil
    % oil eq:
    if vapoil
        oil = (op.pv/dt).*( pvMult.* (bO.* sO  + rv.* bG.* sG) - ...
                              pvMult0.*(bO0.*sO0 + rv0.*bG0.*sG0) ) + ...
                 op.Div(bOvO + rvbGvG);
    else
        oil = (op.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) + op.Div(bOvO);
    end
    if ~isempty(W)
        oil(wc) = oil(wc) - wflux_O;
    end
    eqs{ix} = oil;
    names{ix} = 'oil';
    types{ix} = 'cell';
    ix = ix + 1;
end

if opt.solveForGas
    % gas eq:
    if disgas
        gas = (op.pv/dt).*( pvMult.* (bG.* sG  + rs.* bO.*sO) - ...
                              pvMult0.*(bG0.*sG0 + rs0.*bO0.*sO0 ) ) + ...
                 op.Div(bGvG + rsbOvO);
    else
        gas = (op.pv/dt).*( pvMult.*bG.*sG - pvMult0.*bG0.*sG0 ) + op.Div(bGvG);
    end
    if ~isempty(W)
        gas(wc) = gas(wc) - wflux_G;
    end
    eqs{ix} = gas;
    names{ix} = 'gas';
    types{ix} = 'cell';
    ix = ix + 1;
end

% Polymer eq:
poro = model.rock.poro;
polymer = (op.pv.*(1-fluid.dps)/dt).*(pvMult.*bW.*sW.*c - ...
    pvMult0.*fluid.bW(p0).*sW0.*c0) + (op.pv/dt).* ...
    ( fluid.rhoR.*((1-poro)./poro).*(ads-ads0) ) + op.Div(bWvP);
polymer(wc) = polymer(wc) - wflux_P;

% Fix for (almost) zero water in the well
if isa(polymer, 'ADI')
    is_polymer = strcmpi(primaryVars, 'polymer');
    epsilon = 1.e-8;
    epsilon = sqrt(epsilon)*mean(abs(diag(polymer.jac{is_polymer})));
    bad     = abs(diag(polymer.jac{is_polymer})) < epsilon;
    polymer(bad) = c(bad);
end
bad = value(sW) == 0;
polymer(bad) = c(bad);

names{ix}  = 'polymer';
types{ix} = 'cell';
eqs{ix} = polymer./fluid.cmax; % scale with cmax

% rho = {rhoW, rhoO, rhoG};
% mob = {mobW, mobO, mobG};
% sat = {sW, sO, sG};
% dissolved = model.getDissolutionMatrix(rs, rv);
% eqs = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
%                                      {pW, p, pG}, sat, mob, rho, ...
%                                      dissolved, {c}, ...
%                                      drivingForces);

if ~model.useCNVConvergence
    for i = 1:numel(eqs)
        eqs{i} = eqs{i}.*(dt./op.pv);
    end
end

state.cmax0 = cmax0;  
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

end

%--------------------------------------------------------------------------
% Helper functions
%--------------------------------------------------------------------------

% Effective adsorption, depending of desorption or not
function y = effads(c, cmax, model)
   if model.fluid.adsInx == 2
      y = model.fluid.ads(max(c, cmax));
   else
      y = model.fluid.ads(c);
   end
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

