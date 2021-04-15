function [problem, state] = pressureEquationOilWaterNTPFA(state0, state, model, dt, drivingForces, varargin)

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'staticWells',  false, ...
             'iteration', -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.W;

% assert(isempty(drivingForces.bc) && isempty(drivingForces.src))

s = model.operators;
f = model.fluid;

[p, sW, wellSol] = model.getProps(state, 'pressure', 'water', 'wellsol');
[p0, sW0, wellSol0] = model.getProps(state0, 'pressure', 'water', 'wellsol');
[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

%Initialization of independent variables ----------------------------------

if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [p, wellVars{:}] = initVariablesADI(p, wellVars{:});
    else
        assert(0, 'Backwards solver not supported for splitting');
    end
end
primaryVars = {'pressure', wellVarNames{:}};

% -------------------------------------------------------------------------
sO  = 1 - sW;
sO0 = 1 - sW0;

[krW, krO] = model.evaluateRelPerm({sW, sO});

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

% Modifiy relperm by mobility multiplier (if any)
krW = mobMult.*krW; krO = mobMult.*krO;

% Compute transmissibility
T = s.T.*transMult;

% Gravity contribution
% gdz = model.getGravityGradient();
pW = p;

dp = s.Grad(p);
upc = (dp.value <= 0);
[upcw, upco] = deal(upc);

mobW = krW./f.muW(p);
mobO = krO./f.muO(p);

[Kgrad, model] = model.getPermGradient(p, p0, drivingForces, transMult);

v = Kgrad(p);
vW = -s.faceUpstr(upcw, mobW).*v;
vO = -s.faceUpstr(upco, mobO).*v;


bW = f.bW(p);
bW0 = f.bW(p0);

bO = f.bO(p);
bO0 = f.bO(p0);

rhoW   = bW.*f.rhoWS;
rhoO   = bO.*f.rhoOS;

% These are needed in transport solver, so we output them regardless of
% any flags set in the model.
state = model.storeFluxes(state, vW, vO, []);
state = model.storeUpstreamIndices(state, upcw, upco, []);
if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, []);
    state = model.storeMobilities(state, mobW, mobO, []);
end
% EQUATIONS ---------------------------------------------------------------
% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.
bOvO = s.faceUpstr(upco, bO).*vO;
bWvW = s.faceUpstr(upcw, bW).*vW;


oil = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0) + s.Div(bOvO);
% water:
wat = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + s.Div(bWvW);

if ~isempty(drivingForces.bc) && 0
    bc = drivingForces.bc;
    drivingForces.bc = [];
else
    bc = [];
end

eqTmp = {wat, oil};
[eqTmp, ~, qRes] = addFluxesFromSourcesAndBC(model, eqTmp, ...
                                       {pW, p},...
                                       {rhoW,     rhoO},...
                                       {mobW,     mobO}, ...
                                       {sW, sO}, ...
                                       drivingForces);
wat = eqTmp{1};
oil = eqTmp{2};

if model.outputFluxes
    state = model.storeBoundaryFluxes(state, qRes{1}, qRes{2}, [], drivingForces);
end
[eqs, names, types] = deal({});

% well equations
% Finally, add in and setup well equations
mob = {mobW, mobO};
rho = {rhoW, rhoO};
dissolved = {};
[eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, types, wellSol0, wellSol, wellVars, wellMap, p, mob, rho, dissolved, {}, dt, opt);

eqs{1} = (oil./bO + wat./bW);
names{1} = 'pressure';
types{1} = 'cell';

if ~isempty(bc)
    % Pressure
    isP = strcmpi(bc.type, 'pressure');
    c = sum(model.G.faces.neighbors(bc.face, :), 2);
    if any(isP)
        cP = c(isP);
        eqs{1}(cP) = p(cP) - bc.value(isP);
    end
    % Flux
    if any(~isP)
        cF = c(~isP);
        eqs{1}(cF) = eqs{1}(cF) - bc.value(~isP);
    end
end
eqs{1} = (dt./s.pv).*eqs{1};

state.timestep = dt;
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

for i = 1:numel(W)
    wp = perf2well == i;
    state.wellSol(i).flux = [double(qW(wp)), double(qO(wp))];
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
