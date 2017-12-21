%% Compare MRST results with other simulators
% This example is a two-phase compositional problem where CO2 is injected
% into a mixture of CO2, Methane and Decane. The problem consists of 1000
% cells and is one dimensional for ease of visualization. The problem is
% divided into a large number of timesteps to ensure that the different
% simulators take approximately the same timesteps.
%
% The problem is challenging in terms of fluid physics because the pressure
% is relatively low, which makes the phase behavior highly pressure
% dependent and all components exist in both phases. Since the wells are
% set to bottom hole pressure controls, the fluid volume injected is
% dependent on correctly calculating the mobility and densities in the
% medium.
%
% MRST uses the Peng-Robinson equation of state by default and the Lohrenz,
% Bray and Clark (LBC) correlation to determine viscosities for both
% phases.
mrstModule add compositional deckformat ad-core ad-props
%% Set up model
% MRST includes both natural variables and overall composition. This toggle
% can switch between the modes.
useNatural = true;

pth = getDatasetPath('simplecomp');
fn  = fullfile(pth, 'SIMPLE_COMP.DATA');
% Read deck
deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);
% Set up grid
G = initEclipseGrid(deck);
G = computeGeometry(G);

% Set up rock
rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);
fluid = initDeckADIFluid(deck);
% Define some surface densities
fluid.rhoOS = 800;
fluid.rhoGS = 10;

eos = initDeckEOSModel(deck);
if useNatural
    model = NaturalVariablesCompositionalModel(G, rock, fluid, eos.fluid, 'water', false);
else
    model = OverallCompositionCompositionalModel(G, rock, fluid, eos.fluid, 'water', false);
end
schedule = convertDeckScheduleToMRST(model, deck);

% Manually set the injection composition
[schedule.control.W.components] = deal([0, 1, 0]);
% Injection is pure gas
[schedule.control.W.compi] = deal([1, 0]);

%% Set up initial state
% The problem is defined at 150 degrees celsius with 75 bar initial
% pressure. We set up the initial problem and make a call to the flash
% routines to get correct initial composition.
ncomp = eos.fluid.getNumberOfComponents();

for i = 1:numel(schedule.control.W)
    schedule.control.W(i).lims = [];
end

% Initial conditions
z0 = [0.6, 0.1, 0.3];
T = 150 + 273.15;
p = 75*barsa;

state0 = initCompositionalState(G, p, T, 1, z0, eos);

%% Simulate the schedule
% Note that as the poblem has 500 control steps, this may take some time
% (upwards of 4 minutes).
[ws, states, rep] = simulateScheduleAD(state0, model, schedule);
%% Comparison plots with existing simulators
% The same problem was defined into two other simulators: Eclipse 300
% (which is a commercial simulator) and AD-GPRS (Stanford's research
% simulator). We load in precomputed states from these simulators and
% compare the results.
%
% Note that the shock speed is sensitive to the different tolerances in the
% simulators, which have not been adjusted from the default in either
% simulator. We observe good agreement between all three simulators, with
% the minor differences can likely be accounted for by harmonizing the
% tolerances and sub-timestepping strategy for the different simulators.

lf = get(0, 'DefaultFigurePosition');
h = figure('Position', lf + [0, 0, 350, 0]);
ref = load(fullfile(pth, 'comparison.mat'));
data = {states, ref.statesECL, ref.statesGPRS};
n = min(cellfun(@numel, data));
names = {'MRST', 'E300', 'AD-GPRS'};
markers = {'-', '.', '--'};
cnames = model.EOSModel.fluid.names;

nd = numel(data);
l = cell(nd*ncomp, 1);
for i = 1:nd
    for j = 1:ncomp
        l{(i-1)*ncomp + j} = [names{i}, ' ', cnames{j}];
    end
end
lw = [1, 2, 2];
colors = lines(ncomp);
for step = 1:n
    figure(h); clf; hold on
    for i = 1:numel(data)
        s = data{i}{step};
        comp = s.components;
        if iscell(comp)
            comp = [comp{:}];
        end
        for j = 1:ncomp
            plot(comp(:, j), markers{i}, 'linewidth', lw(i), 'color', colors(j, :));
        end
    end
    legend(l, 'location', 'eastoutside');
    ylim([0, 1]);
    drawnow
end
%% Compare pressure and saturations
% We also plot a more detailed comparison between MRST and E300 to show
% that the prediction of phase behavior is accurate.

colors = lines(ncomp + 2);
for step = 1:n
    figure(h); clf; hold on
    for i = 1:2
        s = data{i}{step};
        if i == 1
            marker = '-';
            linewidth = 1;
        else
            marker = '--';
            linewidth = 2;
        end
        hs = plot(s.s(:, 2), marker, 'color', [0.913, 0.172, 0.047], 'linewidth', linewidth, 'color', colors(1, :));
        p = s.pressure./max(s.pressure);
        hp = plot(p, marker, 'linewidth', linewidth, 'color', colors(2, :));
        comp = s.components;
        if iscell(comp)
            comp = [comp{:}];
        end
        
        if i == 1
            handles = [hs; hp];
        end
        for j = 1:ncomp
            htmp = plot(comp(:, j), marker, 'linewidth', linewidth, 'color', colors(j + 2, :));
            if i == 1
                handles = [handles; htmp];
            end
        end
        if i == 2
            legend(handles, 'sV', 'Normalized pressure', cnames{:}, 'location', 'northoutside', 'orientation', 'horizontal');
        end
    end
    ylim([0, 1]);
    drawnow
end

%% Set up interactive plotting
% Finally we set up interactive plots to make it easy to look at the
% results from the different simulators.

mrstModule add mrst-gui
for i = 1:nd
    figure;
    plotToolbar(G, data{i}, 'plot1d', true);
    title(names{i});
end

%% Plot displacement lines in ternary diagram
figure; hold on
plot([0, 0.5, 1, 0], [0, sqrt(3)/2, 0, 0], 'k')


mapx = @(x, y, z) (1/2)*(2*y + z)./(x + y+ z);
mapy = @(x, y, z) (sqrt(3)/2)*z./(x + y+ z);

colors = parula(numel(states));
for i = 1:5:numel(states)
    C = states{i}.components;
    plot(mapx(C(:, 1), C(:, 2), C(:, 3)), mapy(C(:, 1), C(:, 2), C(:, 3)), '-', 'color', colors(i, :))
    
end
axis off

text(0, 0, model.EOSModel.fluid.names{1}, 'verticalalignment', 'top', 'horizontalalignment', 'right')
text(1, 0, model.EOSModel.fluid.names{2}, 'verticalalignment', 'top', 'horizontalalignment', 'left')
text(0.5, sqrt(3)/2, model.EOSModel.fluid.names{3}, 'verticalalignment', 'bottom', 'horizontalalignment', 'center')


text(mapx(0.5, 0.5, 0), mapy(0.5, 0.5, 0), '0.5', 'verticalalignment', 'top', 'horizontalalignment', 'center')
text(mapx(0, 0.5, 0.5), mapy(0, 0.5, 0.5), '0.5', 'verticalalignment', 'bottom', 'horizontalalignment', 'left')
text(mapx(0.5, 0.0, 0.5), mapy(0.5, 0.0, 0.5), '0.5', 'verticalalignment', 'bottom', 'horizontalalignment', 'right')
