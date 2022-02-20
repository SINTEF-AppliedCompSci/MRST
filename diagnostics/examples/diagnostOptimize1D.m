%% Well Optimization using Adjoints and Flow Diagnostics
% This example demonstrates the diagnostics-based well optimization
% framework on a very simple 1D problem with two injectors and a single
% producer.
%
mrstModule add diagnostics ad-props incomp ad-core

%% Setup model
% We consider a rectangular 10-by-1-by-1 m^3 discretized into 101 cells. We
% divide the reservoir into two regions: One region (x<5) where the
% porosity is set to 0.5 and a region x>5 where the porosity is 1. In the
% center (x = 5), the porosity is set to the average (0.75).
%
% Two injectors are defined. One at the left boundary of the domain and
% another at the right end of the domain. Both are allocated the same
% injection rate, even though they are in different porosity regions.
%
% A single BHP controlled producer is placed in the exact middle of the
% domain, effectively making the reservoir into two distinct regions: The
% low porosity region swept by the first injector and the high porosity
% region swept by the other injector.

% Grid
N = 101; mid = floor(N/2);
G = cartGrid([N, 1, 1], [10, 1, 1]);
G = computeGeometry(G);

% Petrophysics
[lowporo,highporo] = deal(.5,1);
poro = ones(G.cells.num, 1);
poro(1:mid) = lowporo;                 % West part of reservoir
poro(mid+1:end) = highporo;            % East part of reservoir
poro(mid) = (lowporo + highporo)/2 ;   % Center cell
rock = makeRock(G, 1, poro);

T = computeTrans(G, rock);
pv = poreVolume(G, rock);

% Wells
W = verticalWell([], G, rock, 1, 1, [],  ...  % Injector: West
                 'Val', 1*meter^3/day, 'Type', 'rate', ...
                 'Name', 'I1', 'InnerProduct', 'ip_tpf');
             
W = verticalWell(W, G, rock, N, 1, [], ...    % Injector: East
                 'Val', 1*meter^3/day,  'Type', 'rate', ...
                 'Name', 'I2', 'InnerProduct', 'ip_tpf');

W = verticalWell(W, G, rock, mid, 1, [], ...  % Producer: Center
                 'Val', 0, 'Type', 'bhp', ...
                 'Name', 'P', 'InnerProduct', 'ip_tpf');

% Reservoir state
state0 = initResSol(G, 0*barsa, [1 0 0]);
% state0.wellSol = initWellSol(W, 0);

% Reservoir fluid
fluid_ad = initSimpleADIFluid('mu',[1 1 1], 'n', [1 1 1]);

% Set up discrete operators
op = setupOperatorsTPFA(G, rock);
%% Define objective function and plot domain
% The Lorenz coefficient is used here, which will have values between 0
% (homogenous displacement, perfect sweep) and 1 (infinitely hetereogenous
% displacement). Because we want to improve the sweep of the well
% configuration, this is a good choice. The left injector has half the pore
% volume between it and the producer when compared to the other injector,
% and so the optimal configuration should set the ratio between the
% injectors' rates to 1/2.
%
% We at the same time specify the minimum rates per well to a small value,
% which will not be violated at the optimum in this simple  case.
%
% Once the problem definition is complete we plot the wells and porosity on
% the 1D domain.
objective = getObjectiveDiagnostics(G, rock, 'minlorenz');

minRate = W(1).val./1000;
clf
x = G.cells.centroids(:,1);
stairs(x, rock.poro, '-k', 'LineWidth', 2);
grid on
ylim([0, 1.1])
hold on
colors = {'r', 'g', 'b'};
% Plot each well
for i = 1:numel(W)
    c = W(i).cells(1);
    plot(x(c, 1), rock.poro(c), 'O', 'MarkerEdgeColor', 'k', ...
                    'MarkerFaceColor', colors{i}, 'MarkerSize', 8);
end
l = {'Porosity', W.name };
legend(l, 'location', 'southeast')
title('Porosity and well placements')
ylabel('Porosity')
xlabel('X coordinates')

%% Optimize well rates
% We optimize the well rates using a steepest descent-like implementation
% which accounts for well limits and uses adjoints to calculate the
% sensitivites of the objective function. Because the framework uses
% diagnostics for the objective evaluations and adjoints for the
% sensitivities, the cost per function evaluation is low.
%
% We set the tolerance so the iteration has converged when the change
% between iterations is less than 5%. During the optimization, the progress
% will be plotted.
clf
[D_best, W_best, history] = optimizeTOF(G, W, fluid_ad, pv, T, op,...
                                     state0, minRate, objective, ...
                                     'plotProgress', true, ...
                                     'verbose', true, ...
                                     'deltatol', 0.05);

%% Plot the time of flight of the initial and the best well rates
% We plot the time of flight in the domain for each configuration. The
% optimized well configuration shouws equal arrival times from each
% injector, indicating an optimal solution.
clf
pw = @() plotWell(G, W, 'height', 1, 'radius', 0, 'Color', 'red');

subplot(4,1,1:2), title('Time of flight')
x = G.cells.centroids(:,1);
initial = history.D(1).tof(:,1);
optimized = history.D(end).tof(:,1);
hold on
plot(x, [initial, optimized]./max(initial), 'LineWidth', 2)
grid on
legend({'Initial TOF', 'Optimized TOF'});

subplot(4,1,3), title('Initial')
plotCellData(G, history.D(1).tof(:,1),'EdgeColor','k','EdgeAlpha',.1)
[htop, htext] = pw(); %#ok<ASGLU>
set(htext, 'Interpreter', 'None')
axis normal off
view(-4, 40)

subplot(4,1,4), title('Optimized')
colormap winter
plotCellData(G, history.D(end).tof(:,1),'EdgeColor','k','EdgeAlpha',.1)
[htop, htext] = pw();
set(htext, 'Interpreter', 'None')
axis tight off
view(-4, 40)

%% Plot sweep and flow capacity diagrams
% The Lorenz coefficient is really a derived measure from the flow-capacity
% diagram, equal to the integral of the deviation from a perfect, linear
% flow curve where every unit of fluid injected corresponds to a unit of
% recovery. We compute the F and Phi quantities for the base case and the
% optimum and plot both.
%
% At the same time, we plot the sweep diagrams for the same configurations.

[F_0, Phi_0] = computeFandPhi(pv, history.D(1).tof);
[F_end, Phi_end] = computeFandPhi(pv, history.D(end).tof);

clf;
subplot(2,1,1), title('Flow-capacity diagram')
plot([Phi_0, Phi_end], [F_0, F_end], 'linewidth', 2)
legend({'Equal rates', 'Optimized rates'}, 'location', 'SouthEast')
xlabel('\Phi'), ylabel('F'), grid on

subplot(2,1,2), title('Sweep')
[Ev_0, tD_0] = computeSweep(F_0, Phi_0);
[Ev_end, tD_end] = computeSweep(F_end, Phi_end);
plot([tD_0, tD_end], [Ev_0, Ev_end], 'linewidth', 2)
legend({'Equal rates', 'Optimized rates'}, 'location', 'SouthEast')
xlabel('\Phi'), ylabel('F'), grid on, axis tight

%% Copyright notice

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
