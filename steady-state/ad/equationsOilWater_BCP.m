function [problem, state] = equationsOilWater_BCP(state0, state, model, dt, drivingForces, varargin)
% Get linearized problem for oil/water system with black oil-style
% properties

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'iteration', -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.W;
bcp = drivingForces.bcp;
s = model.operators;

% Properties at current timestep
[p, sW, wellSol] = model.getProps(state, 'pressure', 'water', 'wellsol');
% Properties at previous timestep
[p0, sW0, wellSol0] = model.getProps(state0, 'pressure', 'water', 'wellSol');

[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

% Initialize independent variables.
if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [p, sW, wellVars{:}] = initVariablesADI(p, sW, wellVars{:});
    else
        wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
        [p0, sW0, wellVars0{:}] = initVariablesADI(p0, sW0, wellVars0{:}); %#ok
    end
end
% We will solve for pressure, water saturation (oil saturation follows via
% the definition of saturations) and well rates + bhp.
primaryVars = {'pressure', 'sW', wellVarNames{:}};

% Evaluate relative permeability
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
gdz = model.getGravityGradient();

% Evaluate water properties
[vW, bW, mobW, rhoW, pW, upcw] = getFluxAndPropsWater_BO_BCP(model, p, p, sW, krW, T, gdz, bcp);
bW0 = model.fluid.bW(p0);

% Evaluate oil properties
[vO, bO, mobO, rhoO, p, upco] = getFluxAndPropsOil_BO_BCP(model, p, p, sO, krO, T, gdz, [], [], bcp);
bO0 = getbO_BO(model, p0);

if model.outputFluxes
    state = model.storeFluxes(state, vW, vO, []);
end
if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, []);
    state = model.storeMobilities(state, mobW, mobO, []);
    state = model.storeUpstreamIndices(state, upcw, upco, []);
end

% EQUATIONS ---------------------------------------------------------------
% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.
bOvO = s.faceUpstr(upco, bO).*vO;
bWvW = s.faceUpstr(upcw, bW).*vW;

% Conservation of mass for water
water = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + s.Div(bWvW);

% Conservation of mass for oil
oil = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) + s.Div(bOvO);

eqs = {water, oil};
names = {'water', 'oil'};
types = {'cell', 'cell'};

rho = {rhoW, rhoO};
mob = {mobW, mobO};
sat = {sW, sO};


[eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                                 {pW, p}, sat, mob, rho, ...
                                                                 {}, {}, ...
                                                                 drivingForces);
% Finally, add in and setup well equations
[eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, types, wellSol0, wellSol, wellVars, wellMap, p, mob, rho, {}, {}, dt, opt);


% PERIODIC BC & INCOMPRESSIBLE FLUID --------------------------------------

% In the case of periodic boundary conditions AND incompressible fluid, we
% do as in incompTPFA and set a level for the pressure to avoid singular
% Jacobian matrix.
if isempty(W) && ~isempty(bcp)
    % We are dependent on user input. The field 'isIncomp' must be set to
    % true if the fluid is incompressible.
    if isfield(model.fluid, 'isIncomp') && model.fluid.isIncomp
        if eqs{1}.jac{1}(1) > 0
            eqs{1}.jac{1}(1) = 2*eqs{1}.jac{1}(1);
        else
            [~, j] = max(diag(eqs{1}.jac{1}));
            eqs{1}.jac{1}(j,j) = 2*eqs{1}.jac{1}(j,j);
        end
    end
end



problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end
