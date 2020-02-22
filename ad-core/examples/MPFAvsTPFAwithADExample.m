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
if exist('useComp', 'var') && useComp
    % Compositional, two-component
    [f, info] = getCompositionalFluidCase('verysimple');
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
model_mpfa = model.setupStateFunctionGroupings();

mpfa = PhasePotentialDifferenceMPFA(model_mpfa);
model_mpfa.FluxDiscretization = model_mpfa.FluxDiscretization.setStateFunction('PhasePotentialDifference', mpfa);
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
%% Plot the well curves
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
    li = sprintf('%s: Implicit', wn);
    le = sprintf('%s: Explicit', wn);
    % Inconsistent solvers
    qws_te = getWellOutput(wsExplicit, 'qWs', wn);
    qws_ti = getWellOutput(ws, 'qWs', wn);
    figure(h1);
    plot(time/day, -qws_ti, '--', 'linewidth', 2, 'color', c, 'DisplayName', li);
    plot(time/day, -qws_te, '-', 'linewidth', 1, 'color', c, 'DisplayName', le);

    % Consistent solvers
    qws_me = getWellOutput(wsMPFAExplicit, 'qWs', wn);
    qws_mi = getWellOutput(wsMPFA, 'qWs', wn);
    figure(h2);
    plot(time/day, -qws_mi, '--', 'linewidth', 2, 'color', c, 'DisplayName', li);
    plot(time/day, -qws_me, '-', 'linewidth', 1, 'color', c, 'DisplayName', le);
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