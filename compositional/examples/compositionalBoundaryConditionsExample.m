%% Example comparing natural variables and overall composition with boundary conditions
% We simulate injection of CO2 using boundary conditions and compare two
% formulations for the same problem in terms of results and nonlinear
% iterations.
mrstModule add compositional ad-core ad-props
%% Set up problem
nx = 50;
ny = 50;
% Quadratic domain
dims = [nx, ny, 1];
pdims = [1000, 1000, 10];
G = cartGrid(dims, pdims);
G = computeGeometry(G);
% Generate a random Gaussian permeability field
rng(0);
K = logNormLayers(G.cartDims);
K = K(G.cells.indexMap);
rock = makeRock(G, K*milli*darcy, 0.2);

figure; 
plotCellData(G, K);
colorbar;
title('Permeability [mD]');
%% Set up compositional fluid model
% Three-component system of Methane, CO2 and n-Decane
names = {'Methane'  'CarbonDioxide'  'n-Decane'};
% Critical temperatures
Tcrit = [190.5640 304.1282 617.7000];
% Critical pressure
Pcrit = [4599200 7377300 2103000];
% Critical volumes
Vcrit = [9.8628e-05 9.4118e-05 6.0976e-04];
% Acentric factors
acentricFactors = [0.0114 0.2239 0.4884];
% Mass of components (as MRST uses strict SI, this is in kg/mol)
molarMass = [0.0160 0.0440 0.1423];
% Initialize fluid
fluid = CompositionalMixture(names, Tcrit, Pcrit, Vcrit, acentricFactors, molarMass);
%% Set up initial state and boundar conditions
% We inject one pore-volume over 20 years. The right side of the domain is
% fixed at 75 bar.
minP = 75*barsa;
totTime = 20*year;

pv = poreVolume(G, rock);

bc = [];
bc = fluxside(bc, G, 'xmin', sum(pv)/totTime, 'sat', [1, 0]);
bc = pside(bc, G, 'xmax', minP, 'sat', [1, 0]);

bc.components = repmat([0.1, 0.9, 0], numel(bc.face), 1);

flowfluid = initSimpleADIFluid('rho', [1000, 500, 500], ...
                       'mu', [1, 1, 1]*centi*poise, ...
                       'n', [2, 2, 2], ...
                       'c', [1e-5, 0, 0]/barsa);

ncomp = fluid.getNumberOfComponents();
s0 = [1, 0];
state0 = initResSol(G, 1.5*minP, s0);
T = 423.15;
state0.T = repmat(T, G.cells.num, 1);
state0.components = repmat([0.3000 0.1000 0.6000], G.cells.num, 1);

dt = rampupTimesteps(totTime, 75*day, 12);
schedule = simpleSchedule(dt, 'bc', bc);
%% Solve with natural variables
% Natural variables is one possible variable set for compositional
% problems. These variables include saturations, making it easier to limit
% updates for problems with strong relative permeability effects.
model = NaturalVariablesCompositionalModel(G, rock, flowfluid, fluid, 'water', false);

nls = NonLinearSolver('useRelaxation', true);
[~, states, report] = simulateScheduleAD(state0, model, schedule, 'nonlinearsolver', nls);
%% Solve using overall composition
% Another formulation is to use composition variables with no saturations.
% In this formulation, the flash is solved between each nonlinear update.
model_o = OverallCompositionCompositionalModel(G, rock, flowfluid, fluid, 'water', false);

[~, states_o, report_o] = simulateScheduleAD(state0, model_o, schedule, 'nonlinearsolver', nls);
%% Launch interactive plotting
mrstModule add mrst-gui
figure; 
plotToolbar(G, states)
title('Natural variables')

figure; 
plotToolbar(G, states_o)
title('Overall composition')

%% Plot the pressure, gas saturation and CO2 mole fraction for both solvers
% We compare the two results and see that they agree on the solution.
[h1, h2, h3, h4, h5, h6] = deal(nan);
h = figure;
for i = 1:numel(states_o)
    figure(h);
    if i > 1
        delete(h1);
        delete(h2);
        delete(h3);
        delete(h4);
        delete(h5);
        delete(h6);
    end
    subplot(3, 2, 1);
    h1 = plotCellData(G, states{i}.pressure, 'EdgeColor', 'none');
    axis tight; set(gca,'xticklabel',[],'yticklabel',[])
    if i == 1
        title('Natural variables')
        ylabel('Pressure')
    end
    
    subplot(3, 2, 3);
    h3 = plotCellData(G, states{i}.s(:, 2), 'EdgeColor', 'none');
    axis tight; set(gca,'xticklabel',[],'yticklabel',[])
    caxis([0, 1])
    if i == 1
        ylabel('Gas saturation')
    end
    
    subplot(3, 2, 5);
    h5 = plotCellData(G, states{i}.components(:, 2), 'EdgeColor', 'none');
    axis tight; set(gca,'xticklabel',[],'yticklabel',[])
    caxis([0, 1])
    if i == 1
        ylabel(names{2})
    end
    
    subplot(3, 2, 2);
    h2 = plotCellData(G, states_o{i}.pressure, 'EdgeColor', 'none');
    if i == 1
        title('Overall composition')
    end
    axis tight; set(gca,'xticklabel',[],'yticklabel',[])
    
    subplot(3, 2, 4);
    h4 = plotCellData(G, states_o{i}.s(:, 2), 'EdgeColor', 'none');
    axis tight; set(gca,'xticklabel',[],'yticklabel',[])
    caxis([0, 1])
    
    subplot(3, 2, 6);
    h6 = plotCellData(G, states_o{i}.components(:, 2), 'EdgeColor', 'none');
    axis tight off; set(gca,'xticklabel',[],'yticklabel',[])
    caxis([0, 1])
    drawnow
end

%% Plot the nonlinear iteration count
% Different solvers are suitable for different problems. In this case, the
% overall composition model takes slightly fewer iterations.
figure; 
plot([report.Iterations, report_o.Iterations], '-o')
legend('Natural variables', 'Overall composition')
xlabel('Step #')
ylabel('Nonlinear iterations')

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
