%% Combining consistent discretizations with AD-OO
% We follow example 6.1.2 in the MRST book (see
% examples/1ph/showInconsistentTPFA in the book module).
% We create a skewed grid with  two wells where the underlying problem is
% symmetric. An inconsistent discretization of the fluxes may introduce
% asymmetry in the production pattern when injecting a fluid.
mrstModule add ad-core mpfa ad-blackoil compositional ad-props mrst-gui

% Rectangular reservoir with a skew grid.
dims = [41,20];
G = cartGrid(dims,[2,1]);
makeSkew = @(c) c(:,1) + .4*(1-(c(:,1)-1).^2).*(1-c(:,2));
G.nodes.coords(:,1) = 2*makeSkew(G.nodes.coords);
G.nodes.coords(:, 1) = G.nodes.coords(:, 1)*1000;
G.nodes.coords(:, 2) = G.nodes.coords(:, 2)*1000;

G = computeGeometry(G);

% Homogeneous reservoir properties
rock = makeRock(G, 100*milli*darcy, .2);
pv   = sum(poreVolume(G,rock));

% Symmetric well pattern
[ii, jj] = gridLogicalIndices(G);
% Injector + two producers
W = [];
W = addWell(W, G, rock, find(ii == ceil(G.cartDims(1)/2) & jj == G.cartDims(2)), 'comp_i', [1, 0], 'type', 'rate', 'val', pv/year);
W = addWell(W, G, rock, find(ii == G.cartDims(1) & jj == 1), 'comp_i', [1, 0], 'type', 'bhp', 'val', 50*barsa);
W = addWell(W, G, rock, find(ii == 1 & jj == 1), 'comp_i', [1, 0], 'type', 'bhp', 'val', 50*barsa);
%% We can simulate with either immiscible or compositional fluid physics
% The example uses the general simulator framework and as such we can
% easily simulate the same problem with different underlying physics.
gravity reset off;
fluid = initSimpleADIFluid('cR', 1e-8/barsa, 'rho', [1, 1000, 100]);
if ~exist('useComp', 'var')
    useComp = false;
end

if useComp
    % Compositional, two-component
    [f, info] = getBenchmarkMixture('verysimple');
    eos = EquationOfStateModel(G, f);
    model = GenericOverallCompositionModel(G, rock, fluid, eos, 'water', false);
    for i = 1:numel(W)
        W(i).components = info.injection;
    end
    z0 = info.initial;
    state0 = initCompositionalState(G, info.pressure, info.temp, [1, 0], z0, eos);
    W(1).val = 100*W(1).val;
else
    % Immiscible two-phase
    model = GenericBlackOilModel(G, rock, fluid, 'water', true, 'oil', true, 'gas', false);
    state0 = initResSol(G, 1*barsa, [0, 1]);
end
% Schedule
dt = [1; 9; repmat(15, 26, 1)]*day;
schedule = simpleSchedule(dt, 'W', W);
%% Simulate the implicit TPFA base case
[ws, states] = simulateScheduleAD(state0, model, schedule);
%% Simulate implicit MPFA
% The simulator reuses the multipoint transmissibility calculations from
% the MPFA module. We instansiate a special phase potential difference that
% is computed using MPFA instead of the regular two-point difference for
% each face.
mrstModule add mpfa
model_mpfa = setMPFADiscretization(model);
[wsMPFA, statesMPFA] = simulateScheduleAD(state0, model_mpfa, schedule);
%% Simulate explicit MPFA and explicit TPFA
model_exp = setTimeDiscretization(model, 'Explicit', 'initialStep', 0.01*day);
model_mpfa_exp = setTimeDiscretization(model_mpfa, 'Explicit', 'initialStep', 0.01*day);
[wsExplicit, statesExplicit] = simulateScheduleAD(state0, model_exp, schedule);
[wsMPFAExplicit, statesMPFAExplicit] = simulateScheduleAD(state0, model_mpfa_exp, schedule);

%% Plot the results
figure;
plotToolbar(G, states);
title('TPFA')
figure;
plotToolbar(G, statesMPFA)
title('MPFA')
figure;
plotToolbar(G, statesExplicit);
title('TPFA explicit')
figure;
plotToolbar(G, statesMPFAExplicit)
title('MPFA explicit')
%% Plot the producer results
% We note that there are two choices that impact the accuracy of the
% scheme: The choice between a consistent and an inconsistent scheme
% results in bias to the producer that is favorably aligned with the grid
% lines. In addition, the time-stepping and choice of implicitness
% significantly impacts the accuracy of the arrival of the front.
time = cumsum(dt);
h1 = figure; hold on
title('Water production - inconsistent scheme')

h2 = figure; hold on
title('Water production - consistent scheme')
for i = 2:3
    if i == 2
        c = 'r';
    else
        c = 'b';
    end
    wn = W(i).name;
    if useComp
        get = @(ws) -getWellOutput(ws, 'ComponentTotalFlux', wn, 1);
        ti = get(ws);
        te = get(wsExplicit);
        mi = get(wsMPFA);
        me = get(wsMPFAExplicit);
        l = 'CO2 production';
        yl = [0, 0.18];
    else
        rt = {'qWs', 'qOs'};
        qs_te = getWellOutput(wsExplicit, rt, wn);
        qs_ti = getWellOutput(ws, rt, wn);
        qs_me = getWellOutput(wsMPFAExplicit, rt, wn);
        qs_mi = getWellOutput(wsMPFA, rt, wn);

        ti = qs_ti(:, :, 1)./(sum(qs_ti, 3));
        te = qs_te(:, :, 1)./(sum(qs_te, 3));
        mi = qs_mi(:, :, 1)./(sum(qs_mi, 3));
        me = qs_me(:, :, 1)./(sum(qs_me, 3));
        l = 'Water cut';
        yl = [0, 1];
    end
    li = sprintf('%s: Implicit', wn);
    le = sprintf('%s: Explicit', wn);
    % Inconsistent solvers
    figure(h1);
    plot(time/day, ti, '--', 'linewidth', 2, 'color', c, 'DisplayName', li);
    plot(time/day, te, '-', 'linewidth', 1, 'color', c, 'DisplayName', le);
    xlabel('Days simulated');
    ylabel(l);
    ylim(yl);
    % Consistent solvers
    figure(h2);
    plot(time/day, mi, '--', 'linewidth', 2, 'color', c, 'DisplayName', li);
    plot(time/day, me, '-', 'linewidth', 1, 'color', c, 'DisplayName', le);
    xlabel('Days simulated');
    ylabel(l);
    ylim(yl);
end

for i = 1:2
    if i == 1
        figure(h1);
    else
        figure(h2);
    end
    legend('Location', 'northwest');
end
%% Interactive plotting
plotWellSols({ws, wsMPFA, wsExplicit, wsMPFAExplicit}, time, 'datasetnames', {'TPFA implicit', 'MPFA implicit', 'TPFA explicit', 'MPFA explicit'})
%% Plot the front at a chosen time-step
if useComp
    ix = 20;
else
    ix = 10;
end

x = reshape(G.cells.centroids(:, 1), G.cartDims);
y = reshape(G.cells.centroids(:, 2), G.cartDims);

for i = 1:2
    if i == 1
        impl = states;
        expl = statesExplicit;
    else
        impl = statesMPFA;
        expl = statesMPFAExplicit;
    end
    impl = impl{ix};
    expl = expl{ix};
    if useComp
        vi = impl.components(:, 1);
        ei = expl.components(:, 1);
    else
        vi = impl.s(:, 1);
        ei = expl.s(:, 1);
    end
    N        = 20;
    cval     = [.5*movsum(linspace(0,1,N+1),2) 1];
    figure(i); clf; hold on
    colormap(flipud([.7*winter(128).^2+.3; 1 1 1]));
    contourf(x, y, reshape(vi, G.cartDims), cval,'EdgeColor','none');
    contour(x, y, reshape(ei, G.cartDims),  cval(2:end-1), 'k')
    plotGrid(G, 'FaceColor', 'none', 'EdgeColor', 0.8*[1, 1, 1], 'EdgeAlpha', 0.25);

    axis equal tight off;
    for wno = 1:numel(W)
        c = G.cells.centroids(W(wno).cells, :);
        plot(c(1), c(2), 'kO', 'markersize', 8, 'markerFaceColor', 'r')
    end
end

%% Copyright Notice
%
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
