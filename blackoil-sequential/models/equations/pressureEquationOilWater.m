function [problem, state] = pressureEquationOilWater(state0, state, model, dt, drivingForces, varargin)

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'propsPressure', [], ...
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
        [p, wellVars{:}] = model.AutoDiffBackend.initVariablesAD(p, wellVars{:});
    else
        assert(0, 'Backwards solver not supported for splitting');
    end
end
primaryVars = {'pressure', wellVarNames{:}};

p_prop = opt.propsPressure;
otherPropPressure = ~isempty(p_prop);
if ~otherPropPressure
    p_prop = p;
end

% -------------------------------------------------------------------------
sO  = 1 - sW;
sO0 = 1 - sW0;

[krW, krO] = model.evaluateRelPerm({sW, sO});

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p_prop, p0);

% Modifiy relperm by mobility multiplier (if any)
krW = mobMult.*krW; krO = mobMult.*krO;

% Compute transmissibility
T = s.T.*transMult;

% Gravity contribution
gdz = model.getGravityGradient();


% Evaluate water properties
[vW, bW, mobW, rhoW, pW, upcw, dpW] = getFluxAndPropsWater_BO(model, p_prop, sW, krW, T, gdz);
bW0 = f.bW(p0);

% Evaluate oil properties
[vO, bO, mobO, rhoO, pO, upco, dpO] = getFluxAndPropsOil_BO(model, p_prop, sO, krO, T, gdz);
bO0 = getbO_BO(model, p0);

if otherPropPressure
    % We have used a different pressure for property evaluation, undo the
    % effects of this on the fluxes.
    dp_diff = s.Grad(p) - s.Grad(p_prop);
    
    vW = -s.faceUpstr(upcw, mobW).*s.T.*(dpW + dp_diff);
    vO = -s.faceUpstr(upco, mobO).*s.T.*(dpO + dp_diff);
end

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
water = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + s.Div(bWvW);

eqs = {water, oil};

rho = {rhoW, rhoO};
mob = {mobW, mobO};
sat = {sW, sO};

names = {'water', 'oil'};
types = {'cell', 'cell'};
[eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                                 {pW, p}, sat, mob, rho, ...
                                                                 {}, {}, ...
                                                                 drivingForces);
% Finally, add in and setup well equations
dissolved = {};
[eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, types, ...
    wellSol0, wellSol, wellVars, wellMap, p, mob, rho, dissolved, {}, dt, opt);

eqs{1} = (dt./s.pv).*(eqs{1}./bW + eqs{2}./bO);
names{1} = 'pressure';
types{1} = 'cell';
eqs = eqs([1, 3:end]);
names = names([1, 3:end]);
types = types([1, 3:end]);

state.timestep = dt;
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end

%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

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
