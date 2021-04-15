function [problem, state] = transportEquationOilWaterPolymer(state0, state, model, dt, drivingForces, varargin)

opt = struct('Verbose'      , mrstVerbose, ...
             'reverseMode'  , false      , ...
             'scaling'      , []         , ...
             'resOnly'      , false      , ...
             'history'      , []         , ...
             'solveForWater', true       , ...
             'solveForOil'  , false      , ...
             'iteration'    , -1         , ...
             'stepOptions'  , []           );  % Compatibility only

opt = merge_options(opt, varargin{:});
assert(isempty(drivingForces.bc) && isempty(drivingForces.src))

W     = drivingForces.W;
op    = model.operators;
fluid = model.fluid;

% Polymer shear thinning/thickening
usingShear = isfield(fluid, 'plyshearMult');

assert(~(opt.solveForWater && opt.solveForOil));

% Properties at current timestep
[p, sW, c, cmax, wellSol] = model.getProps(state, 'pressure', 'water', ...
    'polymer', 'polymermax', 'wellsol');

% Properties at previous timestep
[p0, sW0, c0, cmax0] = model.getProps(state0, 'pressure', 'water', ...
   'polymer', 'polymermax');

% If timestep has been split relative to pressure, linearly interpolate in
% pressure.
pFlow = p;
if isfield(state, 'timestep')
    dt_frac = dt/state.timestep;
    p = p.*dt_frac + p0.*(1-dt_frac);
end

%Initialization of independent variables ----------------------------------

if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [sW, c] = model.AutoDiffBackend.initVariablesAD(sW, c);
    else
        assert(0, 'Backwards solver not supported for splitting');
    end
end
primaryVars = {'sW', 'polymer'};

clear tmp

% -------------------------------------------------------------------------
sO = 1 - sW;
[krW, krO] = model.evaluateRelPerm({sW, sO});

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

% Modifiy relperm by mobility multiplier (if any)
krW = mobMult.*krW; krO = mobMult.*krO;

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

% Evaluate oil properties
[vO, bO, mobO, rhoO, pO, upco, dpO] = getFluxAndPropsOil_BO(model, p, sO, krO, T, gdz);

gp = op.Grad(p);
Gw = gp - dpW;
Go = gp - dpO;

if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, []);
    state = model.storeMobilities(state, mobW, mobO, mobP);
end

if ~isempty(W)
    
    wflux = sum(vertcat(wellSol.flux), 2);
    perf2well = getPerforationToWellMapping(W);
    wc = vertcat(W.cells);
    
    mobWw = mobW(wc);
    mobOw = mobO(wc);
    
    if usingShear
        % The shear multipliers from the pressure solver are used
        shearMultW = vertcat(state.wellSol.shearMult);
        mobWw = mobWw.*shearMultW;
    end
    
    totMobw = mobWw + mobOw;

    f_w_w = mobWw./totMobw;
    f_o_w = mobOw./totMobw;

    isInj = wflux > 0;
    compWell = vertcat(W.compi);
    compPerf = compWell(perf2well, :);

    f_w_w(isInj) = compPerf(isInj, 1);
    f_o_w(isInj) = compPerf(isInj, 2);

    bWqW = bW(wc).*f_w_w.*wflux;
    bOqO = bO(wc).*f_o_w.*wflux;
    
    % Polymer well equations
    [~, wciPoly, iInxW] = getWellPolymer(W);
    cw        = c(wc);
    cw(iInxW) = wciPoly;
    bWqP      = cw.*bWqW;
    
    % Store well fluxes
    wflux_O = value(bOqO);
    wflux_W = value(bWqW);
    wflux_P = value(bWqP);
    
    for i = 1:numel(W)
        perfind = perf2well == i;
        state.wellSol(i).qOs = sum(wflux_O(perfind));
        state.wellSol(i).qWs = sum(wflux_W(perfind));
        state.wellSol(i).qPs = sum(wflux_P(perfind));
    end

end

% Get total flux from state
flux = sum(state.flux(:,1:2), 2);
vT = flux(model.operators.internalConn);

% Stored upstream indices
[flag_v, flag_g] = getSaturationUpwind(model.upwindType, state, {Gw, Go}, vT, op.T, {mobW, mobO}, op.faceUpstr);
flag = flag_v;


upcw  = flag(:, 1);
upco  = flag(:, 2);

mobOf = op.faceUpstr(upco, mobO);
mobWf = op.faceUpstr(upcw, mobW);
mobPf = op.faceUpstr(upcw, c.*mobP);

if usingShear
    % The shear multipliers from the pressure solver are used
    mobWf = mobWf.*state.shearMult;
    mobPf = mobPf.*state.shearMult;
end

totMob = (mobOf + mobWf);
totMob = max(totMob, sqrt(eps));

% Use concervation equation either for water or for oil
if opt.solveForWater
    f_w = mobWf./totMob;
    bWvW   = op.faceUpstr(upcw, bW).*f_w.*(vT + op.T.*mobOf.*(Gw - Go));

    wat = (op.pv/dt).*(pvMult.*bW.*sW - pvMult0.*fluid.bW(p0).*sW0) + ...
        op.Div(bWvW);
    wat(wc) = wat(wc) - bWqW;
    
    eqs{1} = wat;
    name1  = 'water';
else
    f_o = mobOf./totMob;
    bOvO   = op.faceUpstr(upco, bO).*f_o.*(vT + op.T.*mobWf.*(Go - Gw));
    
    oil = (op.pv/dt).*( pvMult.*bO.*(1-sW) - ...
        pvMult0.*fluid.bO(p0).*(1-sW0) ) + op.Div(bOvO);
    oil(wc) = oil(wc) - bOqO;
    
    eqs{1} = oil;
    name1  = 'oil';
end

% Polymer equations
f_p = mobPf./totMob;
bWvP   = op.faceUpstr(upcw, bW).*f_p.*(vT + op.T.*mobOf.*(Gw - Go));

poro = model.rock.poro;
polymer = (op.pv.*(1-fluid.dps)/dt).*(pvMult.*bW.*sW.*c - ...
    pvMult0.*fluid.bW(p0).*sW0.*c0) + (op.pv/dt).* ...
    ( fluid.rhoR.*((1-poro)./poro).*(ads-ads0) ) + op.Div(bWvP);
polymer(wc) = polymer(wc) - bWqP;

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

names  = {name1, 'polymer'};
types  = {'cell',  'cell'};
eqs{2} = polymer;
% rho = {rhoW, rhoO};
% mob = {mobW, mobO};
% sat = {sW, sO};
% eqs = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
%                                      {pFlow, pFlow}, sat, mob, rho, ...
%                                      {}, {c}, ...
%                                      drivingForces);

eqs{2} = eqs{2}/fluid.cpmax;
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
