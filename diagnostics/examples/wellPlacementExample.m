%% Example demonstrating well placement and optimization
% This example shows the well placement algorithm and/or the well rate
% optimization algorithm on regular cartesian grids. The script optimizes
% the wells and then runs incompressible simulations with different fluid
% properties to evaluate the validity of the optimization.
%
% Depending on the configuration, this example may take some time to run.
% Different options can be changed for different results. The default is to
% optimize the placement of the wells, but not their rates.

mrstModule add deckformat ad-fi diagnostics spe10 mrst-gui

% Start with a five spot-like well setup. Otherwise, place all injectors
% near the middle of the domain.
fiveSpotWells = true;
% Optimize the rates during each step of the well placement
optimizeSubsteps = true;
% Optimize placement. If not, only injection will be optimized
optimizePlacement = true;

% Search radius for movement of wells. Larger radius may converge faster.
wradius = 5;

% Uncomment favored permeability type
% Lognormal
rocktype = 'tarbert';
% Channelized
% rocktype = 'ness';
% Uniform
% rocktype = 'uniform';

% Define grid size. If spe10 data is used, set Nx to a maximum of 60 and Ny
% to a maximum of 220. The default case is a simple 60x60 grid for speed.

Nx = 60;
% Ny = 220;
Ny = 60;
Nz = 1;

switch lower(rocktype)
    case 'ness'
        offset = 65;
        uniformRock = false;
    case 'tarbert'
        offset = 0;
        uniformRock = false;
    case 'uniform'
        offset = 0;
        uniformRock = true;
    otherwise
        error('Unknown setup!')
end

if uniformRock
    % The homogenous test case has optimization of substeps disabled by
    % default
    optimizeSubsteps = false;
    fiveSpotWells = false;
    wradius = 1;
end

%% Set up grid, wells and other required quantities
% We set up wells and other quantities based on the options in the previous
% section. This code is a bit verbose because it allows for different
% options in the same example.

% Grid setup
dims = [Nx, Ny, Nz];
pdims = dims.*[20 10 2].*ft();

G = cartGrid(dims, pdims);
G = computeGeometry(G);
% Rock setup
if uniformRock
    rock.poro = 0.3*ones(G.cells.num, 1);
    rock.perm = 100*milli*darcy*ones(G.cells.num, 1);
else
    rock = SPE10_rock(1:Nx, 1:Ny, (1:Nz) + offset);
    rock.perm = convertFrom(rock.perm, milli*darcy);
    mp = 0.01;
    rock.poro(rock.poro < mp) = mp;
end

pv = poreVolume(G, rock);
T = computeTrans(G, rock);

W = [];
W_prod = [];

% Inject pore volume over 10 years per injector
inRate = sum(pv)/(10*year);
addw =  @(W, i, j, k) verticalWell(W, G, rock, i, j, k, 'Val', inRate, 'Type', 'rate');

midx = ceil(Nx/2);
midy = ceil(Ny/2);

ik = [];
for i = -1:2:1
    for j = -1:2:1
        if fiveSpotWells
            % Set up injector in each corner
            W = addw(W, Nx*(i==1) + 1*(i==-1), Ny*(j==1) + 1*(j==-1), ik);
        else
            % Set up injectors jumbled around the producer
            W = addw(W, midx - 10*i, midy - 10*j, ik);
        end
    end
end
% Producer
W_prod = verticalWell(W_prod, G, rock, midx, midy, [], 'Val', 0, 'Type', 'bhp', 'Name', 'Producer');

nw = numel(W);
% Setup well struct so the first four wells are the injectors, and then the
% targets for the well optimization.
targets = 1:nw;
W = [W; W_prod];

% State
state0 = initResSol(G, 0*barsa, [0 1 0]);

% Fluid, ad system
mu = [1 1 1];
n = [1 1 1];
fluid_ad = initSimpleADIFluid('mu', mu*centi*poise, 'n', n);
sys = initADISystem({'Oil', 'Water'}, G, rock, fluid_ad);

% Minimum fluid injection is 1/4 of the injection rate. This can be
% changed.
minRate = inRate/4;

% Finally, plot everything together
close all
plotCellData(G, log10(rock.perm(:, 1)));
plotWell(G, W);
axis equal tight
%% Set up objective functions and optimize
% The optimization may take some time of both rates and placement is to be
% optimized.

% Lorenz coefficient
objective = getObjectiveDiagnostics(G, rock, 'minlorenz');
% Function handle for pressure solve
solvePressure = @(W, varargin) solveStationaryPressure(G, state0, sys, W, fluid_ad, pv, T, 'objective', objective, varargin{:});

[state, D, grad0] = solvePressure(W);
% Disp initial objective value
grad0.objective.val
optWells = @(W, varargin) optimizeWellPlacementDiagnostics(G, W, rock, objective, targets, D, minRate, state0, fluid_ad, pv, T, sys, varargin{:});
optTOF = @(W, varargin) optimizeTOF(G, W, fluid_ad, pv, T, sys,...
                         state0, minRate, objective, ...
                         'targets', targets, ...
                         'plotProgress', true, varargin{:});

if optimizePlacement
    [W_opt, wellHistory, optHistory] = optWells(W, 'optimizesubsteps', optimizeSubsteps, 'searchRadius', wradius, 'wellSteps', 1, 'plotProgress', true);
    [state_o, D_best, grad] = solvePressure(W_opt, 'linsolve', @mldivide);
else
    [D_best W_best history] = optTOF(W);
    wellHistory = [];
end

%% Plot well paths
% We plot the paths taken by the well placement algorithm. This will only
% plot the porosity if only the injection rates is optimized.

clf
hold on
if uniformRock
    ec = 'k';
else
    ec = 'none';
end
plotCellData(G, rock.poro, 'FaceAlpha', 1, 'EdgeColor', ec, 'edgea', .1)
gc = G.cells.centroids;
colors = {'r', 'g', 'b', 'y'};
for i = 1:numel(targets)
    v = W(targets(i)).cells;
    for j = 1:numel(wellHistory)
        v = [v; wellHistory{j}(i).visited];
    end
    plot3(gc(v, 1), gc(v, 2), -ones(numel(v), 1), ['o--', colors{i}], 'LineWidth', 2, 'MarkerSize', 10, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', colors{i})
end
axis tight equal off
view(0, 89.9)
pc = W_opt(end).cells;
plot3(gc(pc, 1), gc(pc, 2), -ones(numel(v), 1), ['o--', 'w'], 'LineWidth', 2, 'MarkerSize', 10, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'w')
plot3(gc([W_opt.cells], 1), gc([W_opt.cells], 2), -ones(numel([W_opt.cells]), 1), ['X', colors{i}], 'LineWidth', 2, 'MarkerSize', 10, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', colors{i})

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
clear cx
close all
Wells = {W, W_opt};
WellNames = {'Original', 'Optimized placement'};
Phis = [];
Fs = [];
clf
for i = 1:numel(Wells)
    [state, Dx, gradpost] = solvePressure(Wells{i});
    [F, Phi] = computeFandPhi(pv, Dx.tof);
    Lc = computeLorenz(F, Phi);

    subplot(2,1, i)
    plotCellData(G, log10(sum(Dx.tof, 2)))
    if exist('cx', 'var')
        caxis(cx)
    else
        cx = caxis();
    end
    title(sprintf('%s, L_c = %1.2f', WellNames{i}, Lc))
    plotWellsPrint(G, Wells{i}, D)

    axis tight equal off
    colorbar
    Phis = [Phis, Phi];
    Fs = [Fs, F];
end
% Add the baseline optimal case (i.e. a linear function of Phi)
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
legend({fs, 'Improved well placement', 'Idealized displacement'}, 'location', 'South')

%% Run two-phase simulation to validate the results
% To validate that the optimization correlates with the simulated oil
% recovery, we simulate a five year period using two different fluids.
%
% The first fluid contains the same values for both oil and water,
% effectively giving us a red/blue-water problem that should give very good
% oil recovery even for the base case, as the saturation fronts will be
% sharp and the water will displace the "oil" efficiently.
%
% The second fluid object will contain a significant viscosity ratio of 10,
% which makes oil recovery much more challenging because the fronts will be
% smeared. The purpose of this fluid is to validate that the single phase
% optimization procedure still has merit when extended to multiphase
% problems. For simplicity, the problem is incompressible and uses a
% pressure/transport splitting for speed.
%

% Define recovery via pore volume
recov = @(state) 1 - sum(state.s(:,2).*pv)/sum(pv);
close all

% Make a copy of the wells, and explicitly make them two phase wells.
W_initial = W;
W_improved = W_opt;
for i = 1:numel(W_initial)
    W_initial(i).compi = W_initial(i).compi(1:2);
    W_improved(i).compi = W_improved(i).compi(1:2);
end

% Set up fluids
fluids = cell(2,1);
fluids{1} = initSimpleFluid('mu' , [   1,  1]*centi*poise     , ...
                      'rho', [1000, 1000]*kilogram/meter^3, ...
                      'n'  , [   1,   1]);
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
    psolve = @(state, w) incompTPFA(state, G, T, fluid, 'wells', w, 'Verbose', false);

    % Transport solve - implicit transport solver
    tsolve = @(state, w, dT) implicitTransport(state, G, dT, rock, ...
                                                    fluid, 'wells', w, 'Verbose', false);

    % Make convenient function handle for doing both pressure and
    % transport
    solve = @(state, w, dT) tsolve(psolve(state, w), w, dT);

    % Five day timesteps to get a nice animation
    dt = 5*day;
    Time = 0;
    % Five year horizont

    endtime = 5*year;
    while Time < endtime
        % Solve the timestep and store the data
        state_a = solve(state_a, W_initial, dt);
        state_b = solve(state_b, W_improved, dt);

        initial = [initial; state_a]; %#ok
        improved = [improved; state_b]; %#ok
        Time = Time + dt;
        clf
        sa = state_a.s(:,2);
        sb = state_b.s(:,2);

        % Plot the base case
        set(0, 'CurrentFigure', h);
        subplot(2,2,1)
        plotCellData(G, sa);
        title('Oil saturation, base case')
        caxis([0, 1])
        axis tight off

        % Plot the improved case
        subplot(2,2,2)
        plotCellData(G, sb);
        title('Oil saturation, optimized')
        caxis([0, 1])
        axis tight off

        % Plot the historic oil production underway
        subplot(2,2, 3:4)
        plot(cumsum(ones(numel(initial), 1)*dt)/year, ...
            100*arrayfun(recov, [initial, improved]), '-', 'LineWidth', 2)
        xlim([0, endtime/year])
        ylim([0, 100]);
        grid on
        xlabel('Years')
        ylabel('Oil recovered (%)')
        legend('Base case', 'Optimized', 'location', 'southeast')
        title(['Recovery for ', fluid.name, ' fluid']);
        drawnow

        % Ouput to terminal...
        fprintf('Recovery after %s: Initial: %2.1f%%, Modified: %2.1f%%\n',  ...
        formatTimeRange(Time), ...
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
    plotCellData(G, sa)
    plotWellsPrint(G, W, D)
    view(0, 89);
    title(sprintf('Recovery baseline, %s: %2.1f%%', fluids{i}.name, 100*recov(initial(ind))))
    axis tight off
    colorbar
    caxis([0, 1])

    subplot(2, nr, 2 + (i-1)*nr)
    plotCellData(G, sb)
    plotWellsPrint(G, W_opt, D)
    view(0, 89);
    title(sprintf('Recovery improved, %s: %2.1f%%', fluids{i}.name, 100*recov(improved(ind))))
    colorbar
    caxis([0, 1])

    axis tight off
end
%% Plot both problems in the same plot
% We plot the oil recovery for both fluid problems simultanously, showing
% that the optimized placement improves recovery in either case. The effect
% on oil recovery from viscosity differences is also shown to be
% significant.

clf
t = cumsum(ones(numel(states_initial{1}), 1)*dt)/year;
y = [states_initial{1}, states_initial{2}, states_improved{1}, states_improved{2}];
plot(t, 100*arrayfun(recov, y), 'LineWidth', 2)
grid on
xlabel('Years')
ylabel('Oil recovery (%)')
legend({'Five spot, unit fluids', 'Five spot, oil/water', 'Optimized wells, unit fluids', 'Optimized wells, oil/water',}, 'location', 'southeast')
