%% Example demonstrating SPE9 run in Jutul
% This example runs the model from the SPE9 benchmark, which was posed
% twenty years ago to compare contemporary black-oil simulators (Killough,
% 1995). The reservoir is described by a 24 x 25 x 15 grid, having a 10
% degree dipping-angle in the x-direction. The 25 producers initially operate
% at a maximum rate of 1500 STBO/D, which is lowered to 100 STBO/D from day
% 300 to 360, and then raised up again to its initial value until the end
% of simulation at 900 days. The single water injector is set to a maximum
% rate of 5000 STBW/D with a maximum bottom-hole pressure of 4000 psi at
% reference depth. This setup will cause free gas to form after ~100 days
% when the reservoir pressure is reduced below the original saturation
% pressure. The free gas migrates to the top of the reservoir. During the
% simulation most of the wells convert from rate control to pressure
% control. A second problem is a discontinuity in the water-oil capillary
% pressure curve, which may cause difficulties in the Newton solver when
% saturations are changing significantly.
%
% The example simulates the setup using both Jutul and MRST and compares
% the results from these simulators with precomputed simulation results
% from ECLIPSE.  will be slight differences in well responses because
% wells are modelled differently in the two simulations. MRST uses a
% standard Peacemann well model, whereas Jutul uses multisegment wells with
% property tables interpolated and exported from MRST.
%
%   Killough, J. E. 1995. Ninth SPE comparative solution project: A
%   reexamination of black-oil simulation. In SPE Reservoir Simulation
%   Symposium,  12-15 February 1995, San Antonio, Texas. SPE 29110-MS, doi:
%   10.2118/29110-MS

%% Load MRST modules and set up the model based on a standard input deck
mrstModule add ad-core ad-blackoil deckformat jutul

pth = getDatasetPath('spe9');
deck_path  = fullfile(pth, 'BENCH_SPE9.DATA');
problem = initEclipsePackedProblemAD(deck_path);

schedule = problem.SimulatorSetup.schedule;
model = problem.SimulatorSetup.model;

%% Simulate in Jutul
% The simulation does not start automatically. Instead, the routine
% outputs a command which you must paste into an existing Julia session
% that already has JutulDarcy preloaded. Once this simulation is finished,
% you can hit the return button to continue execution in MATLAB. Notice
% that the first time you run the simulation, it will take a long time
% since Julia has to compile necessary code. Once compiled, however, the
% simulation is fast. 
[ws, states] = simulatePackedProblemJutul(problem, 'daemon', false);

%% Simulate in MRST
simulatePackedProblem(problem);
[ws_m, states_m] = getPackedSimulatorOutput(problem);

%% Plot states interactively
mrstModule add mrst-gui
G = model.G;
figure; plotToolbar(G, states, 'field', 'pressure')
view(30, 45), 
plotWell(G, problem.SimulatorSetup.schedule.control(1).W,'FontSize',8);
title('Jutul')

figure; plotToolbar(G, states_m, 'field', 'pressure')
view(30, 45)
plotWell(G, problem.SimulatorSetup.schedule.control(1).W,'FontSize',8);
title('MRST')

%% Plot well responses interactively
n = numel(ws);
T = cumsum(schedule.step.val);

plotWellSols({ws, ws_m}, T, 'datasetnames', {'Jutul', 'MRST'})

%% Compare well responses from ECLIPSE100, Jutul, and MRST
% While ECLIPSE100 and MRST match very well, Jutul gives slightly different
% results because of the different well treatment.
addir = mrstPath('ad-blackoil');
compare = fullfile(addir, 'examples', 'spe9', 'compare');
smry = readEclipseSummaryUnFmt(fullfile(compare, 'SPE9'));

compd = 1:(size(smry.data, 2));
Tcomp =  smry.get(':+:+:+:+', 'YEARS', compd);
T = convertTo(cumsum(schedule.step.val), year);
C = lines(3);

juplot = @(data) plot(T, data, '--', 'linewidth', 3, 'color', C(1, :));
mrstplot = @(data) plot(T, data, '-', 'linewidth', 1, 'color', C(1, :));
compplot = @(data) plot(Tcomp, data, 'ro', 'linewidth', 2, 'color', C(2, :));

figure;
names = {'PROD13', 'PROD18'};
nn = numel(names);
for i = 1:nn

    name = names{i};

    comp = convertFrom(smry.get(name, 'WBHP', compd), psia)';
    mrst = getWellOutput(ws_m, 'bhp', name);
    jutul = getWellOutput(ws, 'bhp', name);

    subplot(nn, 1, i)
    hold on
    mrstplot(mrst);
    compplot(comp);
    juplot(jutul);
    title(name)
    axis tight
    grid on

    xlabel('Time (years)')
    ylabel('Pressure (Pa)')
end
legend({'MRST', 'ECLIPSE', 'Jutul'})
%% Plot all three simulators interactively
[ws_e, T_e] = convertSummaryToWellSols(smry, 'field');
plotWellSols({ws, ws_m, ws_e}, {T, T, T_e}, 'datasetnames', {'Jutul', 'MRST', 'ECLIPSE'})

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
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
