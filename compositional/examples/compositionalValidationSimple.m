%% Validation of MRST against two other simulators
% This example is a two-phase compositional problem in which CO2 is
% injected into a mixture of CO2, Methane, and Decane. The problem consists
% of 1000 cells and is one-dimensional for ease of visualization. The
% problem is divided into a large number of time steps to ensure that the
% different simulators take approximately the same timesteps.
%
% The problem is challenging in terms of fluid physics because the pressure
% is relatively low, which makes the phase behavior highly pressure
% dependent and all components exist in both phases. Since the wells are
% set to bottom-hole pressure controls, the fluid volume injected depends
% on correctly calculating the mobility and densities in the medium.
%
% MRST uses the Peng-Robinson equation of state by default and the Lohrenz,
% Bray, and Clark (LBC) correlation to determine viscosities for both
% phases.
%
% This example is discussed in Section 8.5.1 in the second MRST book:
% Advanced Modelling with the MATLAB Reservoir Simulation Toolbox (MRST),
% Cambridge University Press, 2021.
mrstModule add compositional deckformat ad-core ad-props

%% Set up model
% MRST includes both natural variables and overall composition. This toggle
% can switch between the modes.
if ~exist('useNatural', 'var')
    useNatural = true;
end
[state0, model, schedule, ref] = setupSimpleCompositionalExample(useNatural);
if useNatural
    name = 'Natural';
else
    name = 'Overall'; %#ok<UNRCH>
end
problem = packSimulationProblem(state0, model, schedule, 'simple_comp', 'name', name);

%% Simulate the schedule
% Note that as the problem has 500 control steps, this may take some time
% (upwards of 4 minutes).
simulatePackedProblem(problem);
[ws, states, rep] = getPackedSimulatorOutput(problem);

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
ncomp = model.EOSModel.getNumberOfComponents();

lf = get(0, 'DefaultFigurePosition');
h = figure('Position', lf + [0, 0, 350, 0]);
data = {states, ref.statesECL(2:end), ref.statesGPRS};
n = min(cellfun(@numel, data));
names = {'MRST', 'E300', 'AD-GPRS'};
markers = {'-', '--', '--'};
cnames = model.EOSModel.getComponentNames();

nd = numel(data);
l = cell(nd*ncomp, 1);
for i = 1:nd
    for j = 1:ncomp
        l{(i-1)*ncomp + j} = [names{i}, ' ', cnames{j}];
    end
end
lw = [.5, 2.5, 5];
colors = lines(ncomp);
figure(h)
for step = 1:n % 180 for plot in book
    cla; hold on
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
    legend(l, 'location', 'north', 'numcolumns', 3);
    ylim([0, 1]);
    ylabel('z')
    drawnow
end

%% Compare pressure and saturations
% We also plot a more detailed comparison between MRST and E300 to show
% that the prediction of phase behavior is accurate.

colors = lines(ncomp + 2);
for step = 180 % 1:n
    figure(h); clf; hold on
    for i = 1:2
        s = data{i}{step};
        if i == 1
            marker = '-';
            linewidth = 1;
        else
            marker = '--';
            linewidth = 2.5;
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
                handles = [handles; htmp]; %#ok<AGROW>
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
    plotToolbar(model.G, data{i}, 'plot1d', true);
    title(names{i});
end

%% Plot displacement lines in ternary diagram
figure; hold on
plot([0, 0.5, 1, 0], [0, sqrt(3)/2, 0, 0], 'k')


mapx = @(x, y, z) (1/2)*(2*y + z)./(x + y+ z);
mapy = @(x, y, z) (sqrt(3)/2)*z./(x + y+ z);

colors = parula(numel(states));
for i = 1:20:numel(states)
    C = states{i}.components;
    plot(mapx(C(:, 1), C(:, 2), C(:, 3)), mapy(C(:, 1), C(:, 2), C(:, 3)), '-', 'color', colors(i, :))
end
axis off

text(0, 0, cnames{1}, 'verticalalignment', 'top', 'horizontalalignment', 'right')
text(1, 0, cnames{2}, 'verticalalignment', 'top', 'horizontalalignment', 'left')
text(0.5, sqrt(3)/2, cnames{3}, 'verticalalignment', 'bottom', 'horizontalalignment', 'center')


text(mapx(0.5, 0.5, 0), mapy(0.5, 0.5, 0), '0.5', 'verticalalignment', 'top', 'horizontalalignment', 'center')
text(mapx(0, 0.5, 0.5), mapy(0, 0.5, 0.5), '0.5', 'verticalalignment', 'bottom', 'horizontalalignment', 'left')
text(mapx(0.5, 0.0, 0.5), mapy(0.5, 0.0, 0.5), '0.5', 'verticalalignment', 'bottom', 'horizontalalignment', 'right')

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
