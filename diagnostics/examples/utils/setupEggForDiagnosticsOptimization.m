function problem = setupEggForDiagnosticsOptimization(memberIx)
% setup diagnostics problem for control/trajectory optimization examples

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

assert(memberIx >= 1 && memberIx <=100, 'Realization number must be between 1 and 100, got %d', memberIx);
[G, ~, ~, deck] = setupEGG('realization', memberIx);
G = addBoundingBoxFields(G);
[state0, model, schedule] = initEclipseProblemAD(deck, 'G', G, 'TimestepStrategy', 'none', 'useMex', true);

% make sure gravity is zero
%gboModel = model;
gravity reset off
model.gravity = gravity();

% reset fluid/model to incompressible
model.fluid = convertToIncompFluid(model, 'state', state0);
model.OutputStateFunctions = {};
model.dpMaxRel = inf;

%diagnostics model
pSolver = AGMGSolverAD('tolerance', 1e-8, 'maxIterations', 100, 'keepNumber', G.cells.num);
% phase weights
phw = [0, 1/2.7]; % [no mobile water initially, inverse speed of BL-front]
dmodel = DiagnosticsModel(LinearPressureModel(model), 'state0', state0, 'pSolver', pSolver, 'pressureTol', 1e-8, ...
    'PhaseWeights', phw, 'pvScale', 1.15);

W = schedule.control(1).W;
schedule = simpleSchedule(nan, 'W', W);

% final setup
problem.SimulatorSetup = struct('state0', state0, 'schedule', schedule, 'model', dmodel);
end
