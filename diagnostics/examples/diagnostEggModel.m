%% Ranking of different realizations of the same model
% The Egg model contains a large number of different realizations of the
% same permeability field. We can use flow diagnostics to efficiently rank
% the different models. We recommend that you are familiar with the basic
% concepts of flow diagnostics as demonstreated by the diagnostIntro
% example.
%
% For details on the EggModel and the corresponding ensemble, see
% Jansen, J. D., et al. "The egg model–a geological ensemble for reservoir
% simulation." Geoscience Data Journal 1.2 (2014): 192-195.

%% We first solve the base case using a pressure solver
% Set up pressure solver and solve base case.
mrstModule add deckformat ad-blackoil ad-core ...
   sequential diagnostics example-suite

[G, rock, fluid, deck, state0] = setupEGG('realization', 0);
model = PressureOilWaterModel(G, rock, fluid);
schedule = convertDeckScheduleToMRST(model, deck);
W = schedule.control(1).W;

% Solve base case
state = standaloneSolveAD(state0, model, 1*year, 'W', W);

%% Go through the different cases and store the diagnostics output
% There are 100 additional realizatiosn included, for a total of 101
% different realizations. We compute pressure solutions for all
% realizations and compute diagnostics, including the F-Phi curve.
%
% We also compute the averaged injector and producer tracers, in order to
% assess the uncertainty in sweep and drainage if all the realizations are
% equiprobable.

n = 100;

D_all = cell(n + 1, 1);
rock_all = cell(n + 1, 1);

rock_best = rock;
% Compute base case
D = computeTOFandTracer(state, G, rock, 'wells', W);
rock_all{1} = rock;
D_all{1} = D;
[F, Phi] = deal(zeros(G.cells.num + 1, n+1));
[F(:, 1), Phi(:, 1)] = computeFandPhi(poreVolume(G,rock), D.tof);

ptracer = D.ptracer;
itracer = D.itracer;

t = 0;
hold on
for i = 1:n
    % Compute for all realizations
    fprintf('%d of %d\n', i, n);
    deck = getDeckEGG('realization', i);
    rock  = initEclipseRock(deck);
    rock = compressRock(rock, G.cells.indexMap);
    model = PressureOilWaterModel(G, rock, fluid);
    state = standaloneSolveAD(state0, model, 1*year, 'W', W);

    % Compute diagnostics
    tic();
    D = computeTOFandTracer(state, G, rock, 'wells', W);
    t = t + toc();
    [F(:, i+1), Phi(:, i+1)] = computeFandPhi(poreVolume(G,rock), D.tof);    
    ptracer = ptracer + D.ptracer;
    itracer = itracer + D.itracer;
    rock_all{i+1} = rock;
    D_all{i+1} = D;
end

% Compute averaged injector and producer partitions over all realizations
ptracer = ptracer./(n+1);
itracer = itracer./(n+1);

[val,ppart] = max(ptracer,[],2);
ppart(val==0) = 0;

[val,ipart] = max(itracer,[],2);
ipart(val==0) = 0;
%% Plot the different F-Phi diagrams
% We plot the base case in red, and all others in light gray. We observe
% that the base case is worse than the average, as the deviation from a
% piston-like displacement is large.

figure; hold on
plot(F(:, 2:end), Phi(:, 2:end), 'color', [0.3, 0.3, 0.3])
plot(F(:, 1), Phi(:, 1), 'r', 'linewidth', 2)
%% Iterate over all 
L = zeros(n+1, 1);
for i = 1:n+1
    L(i) = computeLorenz(F(:, i), Phi(:, i));
end

[mv, mix] = min(L);
[Mv, Mix] = max(L);
fprintf('Best case is case #%d with L_c = %1.4f\n', mix, mv);
fprintf('Worst case is case #%d with L_c = %1.4f\n', Mix, Mv);

clf;
subplot(2, 1, 1)
plotCellData(G, rock_all{mix}.perm(:, 1), 'EdgeColor', 'none')
axis tight off
title('Best realization')

subplot(2, 1, 2)
plotCellData(G, rock_all{Mix}.perm(:, 1), 'EdgeColor', 'none')
axis tight off
title('Worst realization')
%% Plot the producer tracer concentrations
% We plot both the base case, and the averaged case where the tracer value
% is equal to the average over all realizations. The uncertain regions that
% change between different wells over the ensamble are colored in gray.
wsign = [W.sign];
figure; 
plotTracerBlend(G, D.ppart, max(D.ptracer, [], 2))
title('Base case production tracers')
view(-50, 70);
axis tight off
plotWell(G, W(wsign > 0), 'fontsize', 0, 'color', 'k')
plotWell(G, W(wsign < 0), 'fontsize', 0, 'color', 'c')

figure;
plotTracerBlend(G, ppart, max(ptracer, [], 2))
title('Average case production tracers');
view(-50, 70);
axis tight off
plotWell(G, W(wsign > 0), 'fontsize', 0, 'color', 'k')
plotWell(G, W(wsign < 0), 'fontsize', 0, 'color', 'c')

%% Plot injector tracers
figure; 
plotTracerBlend(G, D.ipart, max(D.itracer, [], 2))
title('Base case injector tracers');
view(-50, 70);
axis tight off
plotWell(G, W(wsign > 0), 'fontsize', 0, 'color', 'k')
plotWell(G, W(wsign < 0), 'fontsize', 0, 'color', 'c')

figure;
plotTracerBlend(G, ipart, max(itracer, [], 2))
title('Averaged case injector tracers');
view(-50, 70);
axis tight off
plotWell(G, W(wsign > 0), 'fontsize', 0, 'color', 'k')
plotWell(G, W(wsign < 0), 'fontsize', 0, 'color', 'c')

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
