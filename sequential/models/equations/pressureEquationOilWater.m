function [problem, state] = pressureEquationOilWater(state0, state, model, dt, drivingForces, varargin)

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'propsPressure', [], ...
             'staticWells',  false, ...
             'iteration', -1);

opt = merge_options(opt, varargin{:});
s = model.operators;

[p, sW, sO, wellSol] = model.getProps(state, 'pressure', 'water', 'oil', 'wellsol');
[sW0, sO0, wellSol0] = model.getProps(state0, 'water', 'oil', 'wellsol');

[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

%Initialization of independent variables ----------------------------------

if ~opt.resOnly
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode
        [p, wellVars{:}] = model.AutoDiffBackend.initVariablesAD(p, wellVars{:});
    else
        assert(0, 'Reverse mode not supported for splitting solvers');
    end
end
state = model.initStateFunctionContainers(state);
primaryVars = [{'pressure'}, wellVarNames];
p_prop = opt.propsPressure;
if isempty(p_prop)
    state.pressure = p;
else
    state.pressure = p_prop;
end
% -------------------------------------------------------------------------
[b, pv] = model.getProps(state, 'ShrinkageFactors', 'PoreVolume');
[b0, pv0] = model.getProps(state0, 'ShrinkageFactors', 'PoreVolume');

[bW, bO] = deal(b{:});
[bW0, bO0] = deal(b0{:});

[rho, mob, pPhase] = model.getProps(state, 'Density', 'Mobility', 'PhasePressures');

if ~isempty(p_prop)
    state.pressure = p;
end
[phaseFlux, flags] = model.getProps(state, 'PhaseFlux',  'PhaseUpwindFlag');
[vW, vO] = deal(phaseFlux{:});
[upcw, upco] = deal(flags{:});
[mobW, mobO] = deal(mob{:});
[rhoW, rhoO] = deal(rho{:});

% These are needed in transport solver, so we output them regardless of
% any flags set in the model.
state = model.storeFluxes(state, vW, vO, []);
state = model.storeUpstreamIndices(state, upcw, upco, []);
if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, []);
    state = model.storeMobilities(state, mobW, mobO, []);
    state = model.storeDensity(state, rhoW, rhoO, []);
end
% EQUATIONS ---------------------------------------------------------------
% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.
bOvO = s.faceUpstr(upco, bO).*vO;
bWvW = s.faceUpstr(upcw, bW).*vW;

% Oil/water pseudocomponents
oil = (1/dt).*( pv.*bO.*sO - pv.*bO0.*sO0);
water = (1/dt).*( pv.*bW.*sW - pv0.*bW0.*sW0);

eqs = {water, oil};
for i = 1:numel(eqs)
    if isnumeric(eqs{i})
        eqs{i} = model.AutoDiffBackend.convertToAD(eqs{i}, p);
    end
end
sat = {sW, sO};

names = {'water', 'oil'};
types = {'cell', 'cell'};
[eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                                 pPhase, sat, mob, rho, ...
                                                                 {}, {}, ...
                                                                 drivingForces);
% Finally, add in and setup well equations
dissolved = {};
[eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, types, ...
    wellSol0, wellSol, wellVars, wellMap, p, mob, rho, dissolved, {}, dt, opt);

water = s.AccDiv(eqs{1}, bWvW);
oil = s.AccDiv(eqs{2}, bOvO);

eqs{1} = (dt./s.pv).*(water./bW + oil./bO);
names{1} = 'pressure';
types{1} = 'cell';
eqs = eqs([1, 3:end]);
names = names([1, 3:end]);
types = types([1, 3:end]);

state.timestep = dt;
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
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
