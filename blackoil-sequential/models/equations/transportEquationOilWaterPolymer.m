function [problem, state] = transportEquationOilWaterPolymer(state0, ...
    state, model, dt, drivingForces, varargin)

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'scaling', [],...
             'resOnly', false,...
             'history', [],...
             'solveForWater', false, ...
             'solveForOil', true, ...
             'iteration', -1, ...
             'stepOptions', []);  % Compatibility only

opt = merge_options(opt, varargin{:});

W = drivingForces.Wells;
assert(isempty(drivingForces.bc) && isempty(drivingForces.src))

s = model.operators;
f = model.fluid;
G = model.G;

assert(~(opt.solveForWater && opt.solveForOil));

[p, sW, c, cmax, wellSol] = model.getProps(state, 'pressure', 'water', ...
    'polymer', 'polymermax', 'wellsol');

[p0, sW0, c0, cmax0] = model.getProps(state0, 'pressure', 'water', ...
   'polymer', 'polymermax');

wflux = sum(vertcat(wellSol.flux), 2);

%Initialization of independent variables ----------------------------------

if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [sW, c] = initVariablesADI(sW, c);
    else
        assert(0, 'Backwards solver not supported for splitting');
    end
end
primaryVars = {'sW', 'polymer'};

clear tmp

% -------------------------------------------------------------------------
sO = 1 - sW;
[krW, krO] = model.evaluteRelPerm({sW, sO});

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

% Modifiy relperm by mobility multiplier (if any)
krW = mobMult.*krW; krO = mobMult.*krO;

% Compute transmissibility
T = s.T.*transMult;

% Gravity gradient per face
gdz = model.getGravityGradient();

% Evaluate water properties
%[vW, bW, mobW, rhoW, pW, upcw, dpW] = getFluxAndPropsWater_BO(model, p, sW, krW, T, gdz);

% Evaluate water and polymer properties
ads  = effads(c, cmax, model);
ads0 = effads(c0, cmax0, model);
[vW, vP, bW, ~, mobW, mobP, rhoW, pW, upcw, a, dpW] = ...
    getFluxAndPropsWaterPolymer_BO(model, p, sW, c, ads, ...
    krW, T, gdz);

% Evaluate oil properties
[vO, bO, mobO, rhoO, pO, upco, dpO] = getFluxAndPropsOil_BO(model, p, ...
    sO, krO, T, gdz);

gp = s.Grad(p);
Gw = gp - dpW;
Go = gp - dpO;

if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, []);
    state = model.storeMobilities(state, mobW, mobO, mobP);
end

if ~isempty(W)
    perf2well = getPerforationToWellMapping(W);
    wc = vertcat(W.cells);
    
    mobWw = mobW(wc);
    mobOw = mobO(wc);
    mobPw = mobP(wc);
    totMobw = mobWw + mobOw;

    f_w_w = mobWw./totMobw;
    f_o_w = mobOw./totMobw;
    f_p_w = mobPw./totMobw;

    isInj = wflux > 0;
    compWell = vertcat(W.compi);
    compPerf = compWell(perf2well, :);

    f_w_w(isInj) = compPerf(isInj, 1);
    f_o_w(isInj) = compPerf(isInj, 2);
    f_p_w(isInj) = compPerf(isInj, 1);

    bWqW = bW(wc).*f_w_w.*wflux;
    bOqO = bO(wc).*f_o_w.*wflux;
    
    % Polymer injection
    wpolyi = vertcat(W.poly);
    bWqP = bW(wc).*f_p_w.*wpolyi(perf2well).*wflux;

    % Store well fluxes
    wflux_O = double(bOqO);
    wflux_W = double(bWqW);
    wflux_P = double(bWqP);
    
    for i = 1:numel(W)
        perfind = perf2well == i;
        state.wellSol(i).qOs = sum(wflux_O(perfind));
        state.wellSol(i).qWs = sum(wflux_W(perfind));
        state.wellSol(i).qPs = sum(wflux_P(perfind));
    end

end

% Get total flux from state
flux = sum(state.flux, 2);
vT = flux(model.operators.internalConn);

% Stored upstream indices
if model.staticUpwind
    flag = state.upstreamFlag;
else
    flag = multiphaseUpwindIndices({Gw, Go}, vT, s.T, {mobW, mobO}, ...
        s.faceUpstr);
end

upcw  = flag(:, 1);
upco  = flag(:, 2);

mobOf = s.faceUpstr(upco, mobO);
mobWf = s.faceUpstr(upcw, mobW);
mobPf = s.faceUpstr(upcw, mobP);

totMob = (mobOf + mobWf);
totMob = max(totMob, sqrt(eps));

if opt.solveForWater
    f_w = mobWf./totMob;
    bWvW   = s.faceUpstr(upcw, bW).*f_w.*(vT + s.T.*mobOf.*(Gw - Go));
    
    f_p = mobPf./totMob;
    bWvP   = s.faceUpstr(upcw, bW).*f_p.*(vT + s.T.*mobOf.*(Gw - Go));

    wat = (s.pv/dt).*(pvMult.*bW.*sW - pvMult0.*f.bW(p0).*sW0) + ...
        s.Div(bWvW);
    wat(wc) = wat(wc) - bWqW;
    
    poro = model.rock.poro;
    poly = (s.pv.*(1-f.dps)/dt).*(pvMult.*bW.*sW.*c - ...
        pvMult0.*f.bW(p0).*sW0.*c0) + (s.pv/dt).* ...
        ( f.rhoR.*((1-poro)./poro).*(ads-ads0) ) + s.Div(bWvP);
    
    poly(wc) = poly(wc) - bWqP;
    
    eqs{1} = wat;
    eqs{2} = poly;
    names = {'water', 'polymer'};
    types = {'cell', 'cell'};
else
    error('Not implemented for polymer')
    f_o = mobOf./totMob;
    bOvO   = s.faceUpstr(upco, bO).*f_o.*(vT + s.T.*mobWf.*(Go - Gw));

    oil = (s.pv/dt).*( pvMult.*bO.*(1-sW) - pvMult0.*f.bO(p0).*(1-sW0) ) + s.Div(bOvO);
    oil(wc) = oil(wc) - bOqO;
    
    eqs{1} = oil;
    names = {'oil'};
    types = {'cell'};
end
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



