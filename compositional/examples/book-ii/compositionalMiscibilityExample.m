%% Miscible Displacement
% Simulation of miscible processes is a classical application for
% compositional simulators. The displacement efficiency of gas injection
% depends highly on the in situ reservoir pressure and temperature
% conditions. If the displacement front is kept at miscible conditions
% above the minimum miscibility pressure (MMP), the displacement is
% piston-like because there is no surface tension between the residual oil
% and the injected gas. If conditions away from the injection site fall
% below the MMP, the immiscible behavior leads to reduced sweep efficiency,
% as the unfavorable viscosity ratio between the vapor phase formed by
% injection gas and the reservoir oil leads to the formation of viscous
% fingers and lowered recovery.
%
% To illustrate this, we consider injection from the left of a fixed mass
% of CO2 into a 1D reservoir that initially contains CO2, C1, C5, and C12
% in a pure liquid phase found at 150 bars pressure.  We perform seven
% simulations in which the pressure at the producer on the right is varied
% from 70 to 150 bar.
%
% The example is discussed in Section 8.5.4 of the second MRST book:
% Advanced Modelling with the MATLAB Reservoir Simulation Toolbox (MRST),
% Cambridge University Press, 2021. 
mrstModule add compositional ad-core ad-props linearsolvers

%% Notice on Computational Cost
warning('ComputationalCost:High', ...
       ['Please be advised that this example often takes a long time ', ...
        'to run (e.g., more than 2 hours of CPU time)']);
pause(10)

%% Setup the cases as packed problems
pressures = [70, 90, 100, 110, 120, 130, 150];
np        = numel(pressures);
nstep     = 200;
ncell     = 1000;
problems  = cell(1, np);
nls       = NonLinearSolver('useRelaxation', true);
for i = 1:np
    bhp = pressures(i);
    [state0, model, schedule] = ...
        setupMiscibilityTest(bhp*barsa, ncell, nstep, false);
    problems{i} = packSimulationProblem(state0, model, schedule, ...
        'MiscibilityDemo', 'name', sprintf('%d', bhp), 'nonlinearsolver', nls);
end

%% Run the simulations
% By setting sim_multi=true, you can invoke batch processing in the
% background using as many threads as possible. Otherwise, the simulation
% will be performed as normal with output to the command window
if ~exist('sim_multi', 'var')
    sim_multi = false;
end
if sim_multi
    ppm = PackedProblemManager(problems); %#ok<UNRCH>
    ppm.simulateProblemBackground()
    ppm.monitorProgress();
else
    simulatePackedProblem(problems);
end

%% Extract results from the batch simulation
[~, states, reports] = getMultiplePackedSimulatorOutputs(problems);

%% Animation of the displacement front(s)
% For each time step, we plot the gas saturation for all the different
% simulations.
x = model.G.cells.centroids(:, 1);
colors = summer(np);
for step = 1:min(nstep, min(cellfun(@numelData, states)))
    figure(1); clf; hold on
    for i = 1:np
        plot(x, states{i}{step}.s(:, 2), 'color', colors(i, :), 'linewidth', 2);
    end
    drawnow
end

%% Comparison of formulations
% We end the example by comparing the computational efficiency of the
% natural-variables and overall-composition formulation. To accentuate
% differences, we increase the time steps by a factor 5.
bhp = pressures(1);
[state0, modelNat, schedule] = ...
    setupMiscibilityTest(bhp*barsa, ncell, nstep/5, true);
[~, modelOver, ~] = ...
    setupMiscibilityTest(bhp*barsa, ncell, nstep/5, false);

problemNat = packSimulationProblem(state0, modelNat,  schedule, ...
                'MiscibilityDemo', 'nonlinearsolver', nls, ...
                'name', [sprintf('%d', bhp), '_natural_coarsedt']);
problemMole = packSimulationProblem(state0, modelOver, schedule, ...
                'MiscibilityDemo', 'nonlinearsolver', nls, ...
                'name', [sprintf('%d', bhp), '_overall_coarsedt']);
simulatePackedProblem(problemNat);
simulatePackedProblem(problemMole);

%% Extract results and plot the number of iterations
repNat = problemNat.OutputHandlers.reports(:);
repMol = problemMole.OutputHandlers.reports(:);
its = [cellfun(@(x) x.Iterations, repNat)', cellfun(@(x) x.Iterations, repMol)'];
figure;
bar((its),2); legend('Natural variables','Overall composition');


%%
function [state0, model, schedule] = setupMiscibilityTest(bhp, ncell, nstep, useNatural)
    names = {'CO2', 'C1', 'C5', 'C10'};
    Tcrit = [304.7, 190.60, 419.5, 626.0];
    Pcrit = [73.87, 43.04, 37.47, 24.20]*barsa;
    Vcrit = [0.094, 0.098, 0.258, 0.534]/1000;
    acf = [0.225, 0.013, 0.1956, 0.385];
    mw = [44.01, 16.04, 58.12, 134.0]/1000;
    
    bic = zeros(4, 4);
    bic(1, 2) = 0.1;
    bic(4, 2) = 0.041;
    bic = bic + bic';
    cf = CompositionalMixture(names, Tcrit, Pcrit, Vcrit, acf, mw);
    cf = cf.setBinaryInteraction(bic);
    
    G = cartGrid([ncell, 1, 1], [200, 10, 1]);
    G = computeGeometry(G);
    
    rock = makeRock(G, 100*milli*darcy, 0.3);
    
    z0 = [0.33, 0.03, 0.24, 0.40];
    z = [1, 0, 0, 0];
    z = ensureMinimumFraction(z);
    
    T = 373;
    p0 = 150*barsa; % Guesstimate
    cr = 1e-5/barsa;

    rho = [25, 25];
    fluid = initSimpleADIFluid('phases', 'OG', 'cR', cr, 'n', [2, 2], 'rho', rho);
    
    W = [];
    W = addWell(W, G, rock, 1, 'comp_i', [1, 0], 'radius', 0.15, 'type', 'rate', 'val', 50/day);
    W = addWell(W, G, rock, G.cells.num, 'comp_i', [1, 0], 'radius', 0.15, 'type', 'bhp', 'val', bhp);
    time = 2000*day;
    
    for i = 1:numel(W)
        W(i).components = z;
    end

    dt = repmat(time/nstep, nstep, 1);
    if useNatural
        model = GenericNaturalVariablesModel(G, rock, fluid, cf, 'water', false);
    else
        model = GenericOverallCompositionModel(G, rock, fluid, cf, 'water', false);
    end
    schedule = simpleSchedule(dt, 'W', W);
    state0 = initCompositionalState(G, p0, T, [1, 0], z0, model.EOSModel);
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
