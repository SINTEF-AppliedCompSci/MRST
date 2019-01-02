function [problem, state] = transportEquationOilWater(state0, state, model, dt, drivingForces, varargin)

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
W = drivingForces.W;

s = model.operators;
f = model.fluid;

solveAllPhases = opt.solveForWater && opt.solveForOil;

[p, sW, sO, wellSol] = model.getProps(state, 'pressure', 'water', 'oil', 'wellsol');

[p0, sW0, sO0] = model.getProps(state0, 'pressure', 'water', 'oil');

% If timestep has been split relative to pressure, linearly interpolate in
% pressure.
pFlow = p;
if isfield(state, 'timestep')
    dt_frac = dt/state.timestep;
    p = p.*dt_frac + p0.*(1-dt_frac);
end
%Initialization of independent variables ----------------------------------

assert(~opt.reverseMode, 'Backwards solver not supported for splitting');
if solveAllPhases
    if ~opt.resOnly
        [sW, sO] = model.AutoDiffBackend.initVariablesAD(sW, sO);
    end
    primaryVars = {'sW', 'sO'};
    sT = sO + sW;
    [krW, krO] = model.evaluateRelPerm({sW./sT, sO./sT});
else
    if ~opt.resOnly
        sW = model.AutoDiffBackend.initVariablesAD(sW);
    end
    primaryVars = {'sW'};
    sO = 1 - sW;
    sT = ones(size(double(sW)));
    [krW, krO] = model.evaluateRelPerm({sW, sO});
end


clear tmp

% -------------------------------------------------------------------------

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

% Modifiy relperm by mobility multiplier (if any)
krW = mobMult.*krW; krO = mobMult.*krO;

% Compute transmissibility
T = s.T.*transMult;

% Gravity gradient per face
gdz = model.getGravityGradient();

% Evaluate water properties
[vW, bW, mobW, rhoW, pW, upcw, dpW] = getFluxAndPropsWater_BO(model, p, sW, krW, T, gdz);

% Evaluate oil properties
[vO, bO, mobO, rhoO, pO, upco, dpO] = getFluxAndPropsOil_BO(model, p, sO, krO, T, gdz);

gp = s.Grad(p);
Gw = gp - dpW;
Go = gp - dpO;


if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, []);
    state = model.storeMobilities(state, mobW, mobO, []);
    state = model.storeDensity(state, rhoW, rhoO, []);
end

if ~isempty(W)
    wflux = sum(vertcat(wellSol.flux), 2);
    perf2well = getPerforationToWellMapping(W);
    wc = vertcat(W.cells);
    
    mobWw = mobW(wc);
    mobOw = mobO(wc);
    totMobw = mobWw + mobOw;

    f_w_w = sT(wc).*mobWw./totMobw;
    f_o_w = sT(wc).*mobOw./totMobw;

    isInj = wflux > 0;
    compWell = vertcat(W.compi);
    compPerf = compWell(perf2well, :);

    f_w_w(isInj) = compPerf(isInj, 1);
    f_o_w(isInj) = compPerf(isInj, 2);

    bWqW = bW(wc).*f_w_w.*wflux;
    bOqO = bO(wc).*f_o_w.*wflux;

    % Store well fluxes
    wflux_O = double(bOqO);
    wflux_W = double(bWqW);
    
    for i = 1:numel(W)
        perfind = perf2well == i;
        state.wellSol(i).qOs = sum(wflux_O(perfind));
        state.wellSol(i).qWs = sum(wflux_W(perfind));
    end

end

% Get total flux from state
flux = sum(state.flux, 2);
vT = flux(model.operators.internalConn);

rho = {rhoW, rhoO};
mob = {mobW, mobO};
sat = {sW, sO};
G = {Gw, Go};


components = {{sT.*bW,  []}, ...
              {[],      sT.*bO}};
upstr = model.operators.faceUpstr;
[q_phase, q_components] = computeSequentialFluxes(...
    state, G, vT, T, mob, rho, components, upstr, model.upwindType);
[waterFlux, oilFlux] = deal(q_components{:});


if model.outputFluxes
    state = model.storeFluxes(state, q_phase{:}, []);
end

if opt.solveForWater
    wat = (s.pv/dt).*(pvMult.*bW.*sW - pvMult0.*f.bW(p0).*sW0) + s.Div(waterFlux);
    if ~isempty(W)
        wat(wc) = wat(wc) - bWqW;
    end
end

if opt.solveForOil
    oil = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*f.bO(p0).*sO0 ) + s.Div(oilFlux);
    if ~isempty(W)
        oil(wc) = oil(wc) - bOqO;
    end
end

if solveAllPhases
    eqs = {wat, oil};
    names = {'water', 'oil'};    
    types = {'cell', 'cell'};
elseif opt.solveForOil
    eqs = {oil};
    names = {'oil'};
    types = {'cell'};
else
    eqs = {wat};
    names = {'water'};
    types = {'cell'};
end


[eqs, state, src] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                     {pFlow, pFlow}, sat, mob, rho, ...
                                     {}, {}, ...
                                     drivingForces);

if ~model.useCNVConvergence
    for i = 1:(opt.solveForOil+opt.solveForWater)
        eqs{i} = eqs{i}.*(dt./s.pv);
    end
end

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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
