%% Optimizing Well Placement and Rates
% This tutorial shows algorithms for optimized well placement and/or
% optimization of injection rates for regular Cartesian grids. The script
% optimizes the wells and then runs incompressible simulations with
% different fluid properties to evaluate the validity of the optimization.
%
% Depending on the configuration, this example may take some time to run.
% Different options can be changed for different results. The default is to
% optimize the placement of the wells, but not their rates.

mrstModule add ad-core ad-props diagnostics spe10 incomp
close all
gravity reset
%% Parameters that determine the setup
fiveSpotWells = true;          % If false: all injectors near the middle
optimizeSubsteps = true;       % If false: rates not optimized for each step
optimizePlacement = true;      % If false: only iinjection rates optimized
wradius = 5;                   % Search radius for well movement. Larger
                               % values may give faster convergence
rocktype = 'tarbert';          % Values: 'tarbert', 'ness', 'uniform'
[Nx,Ny,Nz] = deal(60,60,1);    % Grid size. SPE: Nx<=60, Ny<=220

switch lower(rocktype)
   case 'ness',
       [offset, uniformRock] = deal(65, false);
   case 'tarbert',
       [offset, uniformRock] = deal(0 , false);
   case 'uniform', 
       [offset, uniformRock] = deal(0 , true );
       [optimizeSubsteps,fiveSpotWells,wradius] = deal(false,false,1);
   otherwise
      error('Unknown setup!')
end

%% Set up reservoir and fluid properties
% This code is a bit verbose because it allows for different options in the
% same example.

% Grid setup
G = cartGrid([Nx, Ny, Nz], [Nx, Ny, Nz].*[20, 10, 2].*ft);
G = computeGeometry(G);

% Rock setup
if uniformRock
    rock.poro = repmat(0.3,             [G.cells.num, 1]);
    rock.perm = repmat(100*milli*darcy, [G.cells.num, 1]);
else
    rock = getSPE10rock(1:Nx, 1:Ny, (1:Nz) + offset);
    rock.poro(rock.poro < 0.01) = 0.01;
end
pv = poreVolume(G, rock);
T = computeTrans(G, rock);

% State

% Fluid
fluid_ad = initSimpleADIFluid('mu', [1,1,1].*centi*poise, 'n', [1,1,1]);

% Build discrete operators
ops = setupOperatorsTPFA(G, rock);


%% Set up wells
% The wells are set to inject one pore volume per injector over a total
% period of ten years. We use the array |targets| to keep track of the
% wells whose rates are to be optimized. If |fiveSpotWells| is true, the
% injectors are set in the corners of the domain. Otherwise, they are
% positioned ten cells away from the producer.
inRate  = sum(pv)/(10*year);                     % Initial rate
minRate = inRate/4;                              % Minimum allowed rate
[midx, midy] = deal(ceil(Nx/2), ceil(Ny/2));     % Midpoint: producer
addw =  @(W, i, j, k) ...
    verticalWell(W, G, rock, i, j, k, 'Val', inRate, 'Type', 'rate', 'refDepth', 0);
W = [];
for i = -1:2:1
   for j = -1:2:1
      if fiveSpotWells
         % Set up injector in each corner
         W = addw(W, Nx*(i==1) + 1*(i==-1), Ny*(j==1) + 1*(j==-1), []);
      else
         % Set up injectors jumbled around the producer
         W = addw(W, midx - 10*i, midy - 10*j, []);
      end
   end
end
targets = 1:numel(W);                   % Indices of wells to be optimized
W = [W; verticalWell([], G, rock, midx, midy, [], ...
                     'Val', 0, 'Type', 'bhp', 'Name', 'P', 'refDepth', 0)];

%% Plot the initial setup
clf
plotCellData(G, log10(rock.perm(:, 1)), 'EdgeAlpha', 0.125);
plotWell(G, W, 'height',.25);
axis tight, set(gca,'DataAspect',[1 1 .01]), view(-30,50)

%% Set up objective functions and optimize
% As objective function for the optimization, we will use the Lorenz
% coefficient, which means that we try to equilibrate the length of all
% flow paths. The optimization may take some time if both rates and
% placement are to be optimized.

% Objective function
objective = getObjectiveDiagnostics(G, rock, 'minlorenz');

% Handle for pressure solver
state0 = initResSol(G, 0*barsa, [0 1 0]);
solvePressure = @(W, varargin) ... 
   solveStationaryPressure(G, state0, ops, W, fluid_ad, pv, T, ...
   'objective', objective, varargin{:});

% Solve pressure and display objective value
[state, D, grad0] = solvePressure(W); grad0.objective.val

% Optimize
if optimizePlacement
    [W_opt, wellHistory, optHistory] = ...
        optimizeWellPlacementDiagnostics(G, W, rock, objective, targets, ...
               D, minRate, state0, fluid_ad, pv,T, ops, ...
               'optimizesubsteps', optimizeSubsteps, ...
               'searchRadius', wradius, 'wellSteps', 1, ...
               'plotProgress', true);
    [state_o, D_best, grad] = ...
        solvePressure(W_opt, 'linsolve', @mldivide);
else
    [D_best, W_best, history] = ...
        optimizeTOF(G, W, fluid_ad, pv, T, ops,state0, minRate, ...
               objective, 'targets', targets, 'plotProgress', true); %#ok<UNRCH>
    wellHistory = [];
end

%% Plot well paths
% We plot the paths taken by the well placement algorithm. This will only
% plot the porosity if only the injection rates is optimized.
clf, hold on
if uniformRock, ec = 'k'; else ec = 'none'; end
plotCellData(G, rock.poro, 'FaceAlpha', 1, 'EdgeColor', ec, 'EdgeAlpha', 0.1)
gc = G.cells.centroids;
colors = {'r', 'w', 'b', 'm'};
for i = 1:numel(targets)
    v = W(targets(i)).cells;
    for j = 1:numel(wellHistory), v = [v; wellHistory{j}(i).visited]; end  %#ok<AGROW>

    plot3(gc(v, 1), gc(v, 2), -ones(numel(v), 1), ...
          ['o--', colors{i}], 'LineWidth', 2, 'MarkerSize', 10, ...
          'MarkerEdgeColor', 'black', 'MarkerFaceColor', colors{i})
end
axis tight equal off
view(0, 89.9)
pc = W_opt(end).cells;
plot3(gc(pc, 1), gc(pc, 2), -ones(numel(v), 1), ['o--', 'w'], ...
      'LineWidth', 2, 'MarkerSize', 10, 'MarkerEdgeColor', 'black', ...
      'MarkerFaceColor', 'w')
plot3(gc([W_opt.cells], 1), gc([W_opt.cells], 2), ...
      -ones(numel([W_opt.cells]), 1), ['X', colors{i}], 'LineWidth', 2, ...
      'MarkerSize', 10, 'MarkerEdgeColor', 'black', ...
      'MarkerFaceColor', colors{i})

%% Plot total travel time for both base and optimized case
% We plot the total travel time (sum of backward and forward time of
% flight) for both the initial and the optimized configuration. To ensure
% that the improvement can be seen, the colorbars are scaled by the base
% case.
%
% If the optimization was successful, the regions with very high total time
% of flight that corresponding to poor sweep should be reduced. The highest
% travel time is colorized dark red.
% 
clear cx, close all
Wells = {W, W_opt};
WellNames = {'Original', 'Optimized placement'};
[Phis, Fs] = deal([]);
clf
for i = 1:numel(Wells)
    [state, Dx, gradpost] = solvePressure(Wells{i});
    [F, Phi] = computeFandPhi(pv, Dx.tof);
    Lc = computeLorenz(F, Phi);

    subplot(2,1, i)
    plotCellData(G, log10(sum(Dx.tof, 2)), 'EdgeColor','none')
    if exist('cx', 'var'), caxis(cx), else cx = caxis();end
    title(sprintf('%s, L_c = %1.2f', WellNames{i}, Lc))
    plotWellsPrint(G, Wells{i}, D)
    axis tight off; colorbar
    
    Phis = [Phis, Phi];                                         %#ok<AGROW>
    Fs = [Fs, F];                                               %#ok<AGROW>
end
% Add the baseline optimal case (i.e., a linear function of Phi)
Phis = [Phis, Phi];
Fs = [Fs, Phi];

%% Plot the storage / flow capacity diagrams
% We also plot the flow capacity / storage capacity diagram (F-Phi) curves
% for both the base and optimized case. As the Lorenz' coefficient is
% defined proportionally to the integral between the idealized curve and
% actual curve, the improvement is easy to see from the resulting plot.
figure;
plot(Phis, Fs, '--', 'LineWidth', 2)
if fiveSpotWells
    fs = 'Five spot';
else
    fs = 'Initial well placement';
end
title('F / Phi-curve before and after optimization');
xlabel('\Phi');
ylabel('F');
grid on
legend({fs, 'Improved well placement', 'Idealized displacement'}, ...
       'location', 'South')

%% Run two-phase simulation to validate the results
% To validate that the optimization correlates with the simulated oil
% recovery, we simulate a five year period using two different fluids.
%
% The first fluid model has idential properties for both phases effectively
% giving us a red/blue-water problem that should give very good oil
% recovery even for the base case, as the saturation fronts will be sharp
% and the water will displace the "oil" efficiently.
%
% The second fluid model has a significant viscosity ratio of 10, which
% makes oil recovery much more challenging because the displacing fluid
% will tend to finger through the resident fluid. The purpose of this setup
% is to validate that the single-phase optimization procedure still has
% merit when extended to multiphase problems. For simplicity, the problem
% is incompressible and uses a pressure/transport splitting for speed.

% Define recovery via pore volume
recov = @(state) 1 - sum(state.s(:,2).*pv)/sum(pv);
close all

% Make a copy of the wells, and explicitly make them two-phase wells.
W_initial = W;
W_improved = W_opt;
for i = 1:numel(W_initial)
    W_initial(i).compi = W_initial(i).compi(1:2);
    W_improved(i).compi = W_improved(i).compi(1:2);
end

% Set up fluids
fluids = cell([2, 1]);
fluids{1} = initSimpleFluid('mu' , [   1,    1]*centi*poise     , ...
                            'rho', [1000, 1000]*kilogram/meter^3, ...
                            'n'  , [   1,    1]);
fluids{1}.name = 'unitary';

fluids{2} = initSimpleFluid('mu' , [   1,  10]*centi*poise     , ...
                            'rho', [1000, 700]*kilogram/meter^3, ...
                            'n'  , [   1,   1]);
fluids{2}.name = 'oil-water';

% Store the different states in cell arrays
states_improved = cell(2,1);
states_initial = cell(2,1);

for i = 1:numel(fluids)
    h = figure;
    fluid = fluids{i};
    state_a = state0;
    state_a.s = state_a.s(:,1:2);
    state_b = state_a;

    initial = [];
    improved = [];

    % Pressure solve - TPFA discretization
    psolve = @(state, w) ...
       incompTPFA(state, G, T, fluid, 'wells', w, 'Verbose', false);

    % Transport solve - implicit transport solver
    tsolve = @(state, w, dT) ...
       implicitTransport(state, G, dT, rock, fluid, ...
                         'wells', w, 'Verbose', false);

    % Make convenient function handle for doing both pressure and
    % transport
    solve = @(state, w, dT) tsolve(psolve(state, w), w, dT);

    % Five day timesteps to get a nice animation
    dt = 5*day;
    Time = 0;
    % Five year horizont

    endtime = 5*year;
    ntschar = numel(formatTimeRange(endtime - dt));
    while Time < endtime
        % Solve the timestep and store the data
        state_a = solve(state_a, W_initial, dt);
        state_b = solve(state_b, W_improved, dt);

        initial = [initial; state_a]; %#ok
        improved = [improved; state_b]; %#ok
        Time = Time + dt;
        sa = state_a.s(:,2);
        sb = state_b.s(:,2);

        % Plot the base case
        set(0, 'CurrentFigure', h);
        subplot(2,2,1), cla
        plotCellData(G, sa, 'EdgeAlpha', 0.125);
        title('Oil saturation, base case')
        caxis([0, 1]), axis tight off

        % Plot the improved case
        subplot(2,2,2), cla
        plotCellData(G, sb, 'EdgeAlpha', 0.125);
        title('Oil saturation, optimized')
        caxis([0, 1]), axis tight off

        % Plot the historic oil production underway
        subplot(2,2, 3:4), cla
        plot(convertTo(cumsum(repmat(dt, [numel(initial), 1])), year), ...
            100*arrayfun(recov, [initial, improved]), '-', 'LineWidth', 2)
        xlim([0, endtime/year]), ylim([0, 100]);
        grid on
        xlabel('Years'), ylabel('Oil recovered (%)')
        legend('Base case', 'Optimized', 'location', 'southeast')
        title(['Recovery for ', fluid.name, ' fluid']);
        drawnow

        % Ouput to terminal...
        fprintf(['Recovery after %-*s: ', ...
                 'Initial: %4.1f%%, Modified: %4.1f%%\n'],  ...
                 ntschar, formatTimeRange(Time), ...
                 100*recov(state_a), ...
                 100*recov(state_b))
    end

    states_initial{i} = initial;
    states_improved{i} = improved;
end

%% Plot the recovery at the half-way point
% If we simulate long enough, almost all oil will certainly be recovered.
% We plot the half-way point after 2.5 years to show the oil saturation
% along with the well positions.
clf
ind = round(2.5*year/dt);

nr = numel(states_initial);
for i = 1:nr
    initial = states_initial{i};
    improved = states_improved{i};
    sa = initial(ind).s(:,2);
    sb = improved(ind).s(:,2);

    subplot(2, nr, 1 + (i-1)*nr)
    plotCellData(G, sa, 'EdgeColor', 'none')
    plotWellsPrint(G, W, D)
    view(0, 89);
    title(sprintf('Recovery baseline, %s: %4.1f%%', ...
                  fluids{i}.name, 100 * recov(initial(ind))))
    axis tight off
    colorbar
    caxis([0, 1])

    subplot(2, nr, 2 + (i-1)*nr)
    plotCellData(G, sb, 'EdgeColor', 'none')
    plotWellsPrint(G, W_opt, D)
    view(0, 89);
    title(sprintf('Recovery improved, %s: %4.1f%%', ...
                  fluids{i}.name, 100 * recov(improved(ind))))
    colorbar
    caxis([0, 1])

    axis tight off
end

%% Plot both problems in the same plot
% We plot the oil recovery for both fluid problems simultanously, showing
% that the optimized placement improves recovery in either case. The effect
% on oil recovery from viscosity differences is also shown to be
% significant.

figure
t = convertTo(cumsum(repmat(dt, [numel(states_initial{1}), 1])), year);
y = [states_initial{1}, states_initial{2}, ...
     states_improved{1}, states_improved{2}];
plot(t, 100*arrayfun(recov, y), 'LineWidth', 2)
grid on
xlabel('Years')
ylabel('Oil recovery (%)')
legend('Five spot, unit fluids',       ...
       'Five spot, oil/water',         ...
       'Optimized wells, unit fluids', ...
       'Optimized wells, oil/water',   ...
       'Location', 'SouthEast')

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
