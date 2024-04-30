%% Combining consistent discretizations with AD-OO
% We follow example 6.1.2 in the MRST book (see
% examples/1ph/showInconsistentTPFA in the book module). We create a
% skewed grid with two wells where the underlying problem is
% symmetric. An inconsistent discretization of the fluxes may introduce
% asymmetry in the production pattern when injecting a fluid.
clear all
close all

mrstModule add ad-core mpfa ad-blackoil compositional ad-props mrst-gui nfvm

% Without twisting the methods should yield the same results
twist = true;
dims = [41, 20];
G = cartGrid(dims, [2, 1]);
if twist
    makeSkew = @(c) c(:, 1) + .4 * (1 - (c(:, 1) - 1).^2) .* (1 - c(:, 2));
    G.nodes.coords(:, 1) = 2 * makeSkew(G.nodes.coords);
    G.nodes.coords(:, 1) = G.nodes.coords(:, 1) * 1000;
    G.nodes.coords(:, 2) = G.nodes.coords(:, 2) * 1000;
    G = twister(G);
end
G = computeGeometry(G);

% Homogeneous reservoir properties
rock = makeRock(G, 100*milli*darcy, .2);
pv = sum(poreVolume(G, rock));

% Symmetric well pattern
if any(contains(G.type, 'cartGrid'))
    [ii, jj] = gridLogicalIndices(G);
    c1 = find(ii == ceil(G.cartDims(1) / 2) & jj == G.cartDims(2));
    c2 = find(ii == G.cartDims(1) & jj == 1);
    c3 = find(ii == 1 & jj == 1);
else
    offset = 1e-2;
    xmin = min(G.nodes.coords(:, 1:2));
    xmax = max(G.nodes.coords(:, 1:2));
    c1 = findEnclosingCell(G, [0.5 * (xmin(1) + xmax(1)), xmax(2) - offset]);
    c2 = findEnclosingCell(G, xmin+offset);
    c3 = findEnclosingCell(G, [xmax(1) - offset, xmin(2) + offset]);
end
% Injector + two producers
r = 0.005;
W = [];
W = addWell(W, G, rock, c1, 'comp_i', [1, 0], 'type', 'rate', 'val', pv/year, 'radius', r);
W = addWell(W, G, rock, c2, 'comp_i', [1, 0], 'type', 'bhp', 'val', 50*barsa, 'radius', r);
W = addWell(W, G, rock, c3, 'comp_i', [1, 0], 'type', 'bhp', 'val', 50*barsa, 'radius', r);

%% We can simulate with either immiscible or compositional fluid physics
% The example uses the general simulator framework and as such we can
% easily simulate the same problem with different underlying physics.

gravity reset off;

fluid = initSimpleADIFluid('cR', 1e-8/barsa, 'rho', [1, 1000, 100]);
useComp = true;

if useComp
    % Compositional, two-component
    [f, info] = getCompositionalFluidCase('verysimple');
    eos = EquationOfStateModel(G, f);
    model = GenericOverallCompositionModel(G, rock, fluid, eos, 'water', false);
    for i = 1:numel(W)
        W(i).components = info.injection;
    end
    z0 = info.initial;
    state0 = initCompositionalState(G, info.pressure, info.temp, [1, 0], z0, eos);
    W(1).val = 100 * W(1).val;
else
    % Immiscible two-phase
    model = GenericBlackOilModel(G, rock, fluid, 'water', true, 'oil', true, 'gas', false);
    state0 = initResSol(G, 1*barsa, [0, 1]);
end
% Schedule
dt = [1; 9; repmat(15, 26, 1)] * day;
schedule = simpleSchedule(dt, 'W', W);

%% Simulate the implicit TPFA base case
disp('TPFA implicit')
[wsTPFA, statesTPFA] = simulateScheduleAD(state0, model, schedule);
plotFinalPressure(G, statesTPFA, 'TPFA')

%% Simulate implicit AvgMPFA
disp('AvgMPFA implicit')
ratio = [];
model_avgmpfa = setAvgMPFADiscretization(model, 'myRatio', ratio);
[wsAvgMPFA, statesAvgMPFA] = simulateScheduleAD(state0, model_avgmpfa, schedule);
plotFinalPressure(G, statesAvgMPFA, 'AvgMPFA')

%% Simulate implicit NTPFA
disp('NTPFA implicit')
ratio = [];
model_ntpfa = setNTPFADiscretization(model, 'myRatio', ratio);
[wsNTPFA, statesNTPFA] = simulateScheduleAD(state0, model_ntpfa, schedule);
plotFinalPressure(G, statesNTPFA, 'NTPFA')

%% Simulate implicit MPFA
% The simulator reuses the multipoint transmissibility calculations from
% the MPFA module. We instantiate a special phase potential difference that
% is computed using MPFA instead of the regular two-point difference for
% each face.
disp('MPFA implicit')
model_mpfa = setMPFADiscretization(model);
[wsMPFA, statesMPFA] = simulateScheduleAD(state0, model_mpfa, schedule);
plotFinalPressure(G, statesMPFA, 'MPFA')

return

%% Collect implicit solver data
wsImplicit{1} = wsTPFA;
wsImplicit{2} = wsAvgMPFA;
wsImplicit{3} = wsNTPFA;
wsImplicit{4} = wsMPFA;
statesImplicit{1} = statesTPFA;
statesImplicit{2} = statesAvgMPFA;
statesImplicit{3} = statesNTPFA;
statesImplicit{4} = statesMPFA;
names = {'TPFA', 'AvgMPFA', 'NTPFA', 'MPFA'};

%% Simulate using explicit time discretization
models_exp{1} = setTimeDiscretization(model, 'Explicit');
models_exp{2} = setTimeDiscretization(model_avgmpfa, 'Explicit');
models_exp{3} = setTimeDiscretization(model_ntpfa, 'Explicit');
models_exp{4} = setTimeDiscretization(model_mpfa, 'Explicit');

for k = 1:numel(models_exp)
    str = sprintf('%s explicit', names{k});
    fprintf('%s\n', str);
    [wsExplicit{k}, statesExplicit{k}] = simulateScheduleAD(state0, models_exp{k}, schedule); %#ok
    plotFinalPressure(G, statesExplicit{k}, str);
end

%% Plot the producer results
% We note that there are two choices that impact the accuracy of the
% scheme: The choice between a consistent and an inconsistent scheme
% results in bias to the producer that is favorably aligned with the grid
% lines. In addition, the time-stepping and choice of implicitness
% significantly impacts the accuracy of the arrival of the front.
time = cumsum(dt);
color = lines(2*numel(wsImplicit));

h1 = figure; hold on; grid on
title('Water production inconsistent scheme well');

h2 = figure; hold on; grid on
title('Water production consistent schemes well');

for iwell = 2:3
    wn = W(iwell).name;

    if useComp
        get = @(ws) -getWellOutput(ws, 'ComponentTotalFlux', wn, 1);
        for k = 1:numel(wsExplicit)
            out_i{k} = get(wsImplicit{k}); %#ok
            out_e{k} = get(wsExplicit{k}); %#ok
        end
        l = 'CO2 production';
        yl = [0, 0.18];
    else
        rt = {'qWs', 'qOs'};
        get = @(ws) getWellOutput(ws, rt, wn);
        for k = 1:numel(wsExplicit)
            qs_e = get(wsExplicit{k});
            qs_i = get(wsImplicit{k});
            out_i{k} = qs_i(:,:,1) ./ (sum(qs_i, 3)); %#ok
            out_e{k} = qs_e(:,:,1) ./ (sum(qs_e, 3)); %#ok
        end
        l = 'Water cut';
        yl = [0, 1];
    end

    iname = @(name) sprintf('%s %s', name, wn);
    ename = @(name) sprintf('%s explicit %s', name, wn);

    % Inconsistent solvers
    figure(h1);
    plot(time/day, out_i{1}, '-', 'linewidth', 2, 'color', color(iwell-1,:), 'DisplayName', iname('TPFA'));
    plot(time/day, out_e{1}, '--', 'linewidth', 1, 'color', color(iwell-1,:), 'DisplayName', ename('TPFA'));
    xlabel('Days simulated');
    ylabel(l);
    ylim(yl);

    % Consistent solvers
    figure(h2);
    for k = 2:numel(wsExplicit)
        cidx = k-1 + (iwell-2)*numel(wsExplicit);
        plot(time/day, out_i{k}, '-', 'linewidth', 2, 'color', color(cidx,:), 'DisplayName', iname(names{k}));
        plot(time/day, out_e{k}, '--', 'linewidth', 1, 'color', color(cidx,:), 'DisplayName', ename(names{k}));
    end
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
namesExplicit = cellfun(@(name) sprintf('%s explicit', name), names, 'un', false);
plotWellSols({wsImplicit{:}, wsExplicit{:}}, time, ...
             'datasetnames', {names{:}, namesExplicit{:}}); %#ok

%% Plot the front at a chosen time-step
if useComp
    timestep = 20;
else
    timestep = 10;
end

x = reshape(G.cells.centroids(:, 1), G.cartDims);
y = reshape(G.cells.centroids(:, 2), G.cartDims);

for k = 1:numel(statesImplicit)
    impl = statesImplicit{k}{timestep};
    expl = statesExplicit{k}{timestep};

    if useComp
        vi = impl.components(:, 1);
        ei = expl.components(:, 1);
    else
        vi = impl.s(:, 1);
        ei = expl.s(:, 1);
    end

    N = 20;
    cval = [.5 * movsum(linspace(0, 1, N + 1), 2), 1];
    figure(k);
    clf;
    hold on
    colormap(flipud([.7 * winter(128).^2 + .3; 1, 1, 1]));
    contourf(x, y, reshape(vi, G.cartDims), cval, 'EdgeColor', 'none');
    contour(x, y, reshape(ei, G.cartDims), cval(2:end - 1), 'k')
    plotGrid(G, 'FaceColor', 'none', 'EdgeColor', 0.8*[1, 1, 1], 'EdgeAlpha', 0.25);

    axis equal tight off;
    for wno = 1:numel(W)
        c = G.cells.centroids(W(wno).cells, :);
        plot(c(1), c(2), 'kO', 'markersize', 8, 'markerFaceColor', 'r')
    end
    title(names{k});
end


%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.
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
