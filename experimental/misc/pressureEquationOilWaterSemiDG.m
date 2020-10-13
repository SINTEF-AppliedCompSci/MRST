function [problem, state] = pressureEquationOilWaterSemiDG(state0, state, model, dt, drivingForces, varargin)

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'propsPressure', [], ...
             'staticWells',  false, ...
             'iteration', -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.W;

% assert(isempty(drivingForces.bc) && isempty(drivingForces.src))

op = model.operators;
f  = model.fluid;
disc = model.disc;
G  = model.G;

[p, sWdof, sOdof, wellSol] = model.getProps(state, 'pressure', 'water', 'oil', 'wellsol');
[p0, sW0dof, sO0dof, wellSol0] = model.getProps(state0, 'pressure', 'water', 'oil', 'wellsol');

[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

%Initialization of independent variables ----------------------------------

if ~opt.resOnly
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode
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

sW = state.s(:,1);
sO = state.s(:,2);
sW0 = state0.s(:,1);
sO0 = state0.s(:,2);
[krW, krO] = model.evaluateRelPerm({sW, sO});

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p_prop, p0);
pvMult  = expandSingleValue(pvMult , G);
pvMult0 = expandSingleValue(pvMult0, G);
mobMult = expandSingleValue(mobMult, G);

% Modifiy relperm by mobility multiplier (if any)
krW = mobMult.*krW; krO = mobMult.*krO;

% Compute transmissibility
T = op.T.*transMult;

% Gravity contribution
gdz = model.getGravityGradient();


% Evaluate water properties
[vW, bW, mobW, rhoW, pW, upcw, dpW, muW] = getPropsWater_DG(model, p_prop, T, gdz, mobMult);
bW0 = f.bW(p0);

% Evaluate oil properties
[vO, bO, mobO, rhoO, pO, upco, dpO, muO] = getPropsOil_DG(model, p_prop, T, gdz, mobMult);
bO0 = getbO_BO(model, p0);

faces = find(disc.internalConn);
x = G.faces.centroids(faces,:);

upCellsW = G.faces.neighbors(faces,1).*upcw + G.faces.neighbors(faces,2).*(~upcw); 
sWf = disc.evaluateDGVariable(x, upCellsW, state, sWdof);

vW = vW(sWf, 1, upCellsW);

upCellsO = G.faces.neighbors(faces,1).*upco + G.faces.neighbors(faces,2).*(~upco); 
sOf = disc.evaluateDGVariable(x, upCellsO, state, sOdof);
vO = vO(sOf, 1, upCellsO);

if otherPropPressure
    % We have used a different pressure for property evaluation, undo the
    % effects of this on the fluxes.
    dp_diff = op.Grad(p) - op.Grad(p_prop);
    vW = -mobW(sWf, 1, upCellsW).*op.T.*(dpW + dp_diff);
    vO = -mobO(sOf, 1, upCellsO).*op.T.*(dpO + dp_diff);
end

% These are needed in transport solver, so we output them regardless of
% any flags set in the model.
state = model.storeFluxes(state, vW, vO, []);
state = model.storeUpstreamIndices(state, upcw, upco, []);
mobW = mobW(sW, 1, (1:G.cells.num)');
mobO = mobO(sO, 1, (1:G.cells.num)');

if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, []);
    state = model.storeMobilities(state, mobW, mobO, []);
end
% EQUATIONS ---------------------------------------------------------------
% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.
bOvO = op.faceUpstr(upco, bO).*vO;
bWvW = op.faceUpstr(upcw, bW).*vW;


oil = (op.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0) + op.Div(bOvO);
% water:
water = (op.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + op.Div(bWvW);

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

eqs{1} = (dt./op.pv).*(eqs{1}./bW + eqs{2}./bO);
names{1} = 'pressure';
types{1} = 'cell';
eqs = eqs([1, 3:end]);
names = names([1, 3:end]);
types = types([1, 3:end]);

state.timestep = dt;
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end

function v = expandSingleValue(v,G)

    if numel(double(v)) == 1
        v = v*ones(G.cells.num,1);
    end
    
end

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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
