%% Numerical Accuracy for Different Discretizations
% This example is a continuation of the compositionalValidationSimple
% example, in which we compare the numerical accuracy of simulations with
% explicit and implicit discretizations for the same case. This example is
% discussed in Section 8.5.2 of the second MRST book: Advanced Modelling
% with the MATLAB Reservoir Simulation Toolbox (MRST), Cambridge University
% Press, 2021.
mrstModule add ad-core compositional deckformat ad-props

%% Notice on Computational Cost
warning('ComputationalCost:High', ...
       ['Please be advised that this example often takes a long time ', ...
        'to run (e.g., more than 1 hour of CPU time)']);
pause(10)

%% Setup simulation models
% We set up simulation models for both the natural-variable and overall-
% composition formulation and organize these as so-called packed simulation
% cases so that results are stored to disk, once computed. This enables
% easy restart, should the simulations be interrupted. Once the simulations
% are fully completed and stored on disk, you do not have to repeat them if
% you decide to run the script again to, e.g., look at the results.
[state0, modelOverall, schedule] = setupSimpleCompositionalExample(false);
[~,      modelNatural, ~       ] = setupSimpleCompositionalExample(true);

packer = @(model, name, varargin) ...
    packSimulationProblem(state0, model, schedule, 'compAccExample', ...
                            'name', name, varargin{:});
natural = packer(modelNatural, 'Natural');
overall = packer(modelOverall, 'Overall');

%% Add an extra explicit simulation for the overall compositional case
modelExpl = setTimeDiscretization(modelOverall, 'explicit');
explicit = packer(modelExpl, 'Explicit');

%% Simulate the problems
% Notice that these simulation will take long time to run, primarily
% because of the time-step restriction in the explicit scheme, which forces
% the simulator to take a large number of small time steps.
problems = {explicit, natural, overall};
np = numel(problems);
simulatePackedProblem(problems);

%% Load the computed data
[~, states, reports, names] = getMultiplePackedSimulatorOutputs(problems);

%%
ncomp = modelExpl.EOSModel.getNumberOfComponents();

lf = get(0, 'DefaultFigurePosition');
h = figure('Position', lf + [0, 0, 350, 0]);
n = min(cellfun(@numel, states));
markers = {'-', '--', '--'};
cnames = modelExpl.EOSModel.getComponentNames();

nd = numel(states);
l = cell(nd*ncomp, 1);
for i = 1:nd
    for j = 1:ncomp
        l{(i-1)*ncomp + j} = [names{i}, ' ', cnames{j}];
    end
end
lw = [1, 2.5, 5];
colors = lines(ncomp);
figure(h); clf
for step = 1:10:n % set n=110 for plot in the MRST book
    cla; hold on
    for i = 1:numel(states)
        s = states{i}{step};
        comp = s.components;
        for j = 1:ncomp
            plot(comp(:, j), markers{i}, 'linewidth', lw(i), 'color', colors(j, :));
        end
    end
    legend(l, 'location', 'north', 'numcolumns', 3);
    ylim([0, 1]);
    ylabel('z')
    drawnow
end
set(gca,'FontSize',12)

%% Plot gas saturation
colors = lines(numel(states));
lwi = 1.5;
figure(h); clf
for step = 1:10:n % set n=110 for plot in the MRST book
    cla; hold on
    for i = 1:numel(states)
        s = states{i}{step};
        sg = s.s(:, 2);
        if i == 3
            st = '--'; lwi = 4;
        else
            st = '-';
        end
        plot(sg, st, 'linewidth', lwi, 'color', colors(i, :));
    end
    legend(names, 'location', 'north', 'numcolumns', 3);
    ylim([0, 1]);
    ylabel('z')
    drawnow
end

%% Plot properties: viscosity, density, gas saturation
model = modelOverall.validateModel();
lwi = 1.5;
colors = lines(numel(states));
figure(h); clf
for step = 1:10:n % set n=180 for plot in the MRST book
    clf; hold on
    for i = 1:2
        s = states{1}{step};
        mask = s.L == (1 - double(i == 1));
        mu = value(model.getProp(s, 'Viscosity'));
        mu = mu(:, i)/(centi*poise);
        mu(mask) = nan;

        rho = value(model.getProp(s, 'Density'));
        rho = rho(:, i);
        rho(mask) = nan;

        yyaxis left
        plot(mu, 'linewidth', lwi);
        ylabel('Viscosity [cP]')

        yyaxis right
        plot(rho, 'linewidth', lwi);
        ylabel('Mass density [kg/m^3]');
        
    end
    yyaxis right
    yl = ylim();
    plot(yl(2)*s.s(:, 2), '-k', 'linewidth', lwi);
    plot(yl(2)*s.components(:, 2), '-r', 'linewidth', lwi);
    legend('\mu_l', '\mu_v', '\rho_l', '\rho_v', 'S_g', 'z (CO_2)', ...
            'location', 'east', 'numcolumns', 3);
    drawnow
end

%% Plot the cumulative number of iterations for the implicit schemes
natIts = cellfun(@(x) x.Iterations, reports{2});
molIts = cellfun(@(x) x.Iterations, reports{3});
figure;
plot(cumsum([natIts, molIts]))
legend(names(2:end))
title('Cumulative number of iterations')

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
