%% More Accurate Gas Production by Modeling Separators
% This basic example demonstrates how you can set up the simulation to
% include the modelling of one or more separator units as part of the
% production facilities. Each separator is set to operate at a specific
% pressure and temperatur condition and takes the input stream of fluids,
% flashes it to determine the new composition at vapor-liquid equilibrium,
% and direct the new vapor and liquid content into separate output streams.
% In case of multiple separators, these are typically specified as a tree
% structure.
%
% The example is discussed in Section 8.5.3 of the second MRST book:
% Advanced Modelling with the MATLAB Reservoir Simulation Toolbox (MRST),
% Cambridge University Press, 2021. 
mrstModule add ad-core compositional deckformat ad-props

%% Setup the base case
% The base case consists of a homogeneous reservoir described on a 11x11x10
% Cartesian grid. The fluid system consists of six hydrocarbon components
% (C1,C3,C6,C10,C15,C20) as described in the 5th SPE Comparative Solution
% Project. The reservoir is produced from a single injector placed in the
% middle of the domain, which operates a fixed rate that extracts a
% significant volume so that a vapor phase eventually is formed as the
% pressure drops.
%
% By default, MRST uses a simple rule for produced liquids: The surface
% density specified in the fluid object (e.g., fluid.rhoGS for the default
% vapor phase) is used to set liquid and vapor rates. Components are
% separated into liquid and vapor streams according to their phase
% fractions at the standard conditions specified by the facility model.

% Reservoir geometry and petrophysical properties
G    = cartGrid([11, 11, 10], [100, 100, 30]);
G    = computeGeometry(G);
rock = makeRock(G, 0.1*darcy, 0.3);

% Fluid system: SPE5 mixture combined with standard black-oil 
[mix, info] = getBenchmarkMixture('spe5');
fluid = initSimpleADIFluid('blackoil', false, 'phases', 'OG', ...
                            'n', [2, 2], 'rho', [100, 1]);
model = GenericOverallCompositionModel(G, rock, fluid, mix, 'water', false);

% Single well in the middle of the reservoir
W = verticalWell([], G, rock, 6, 6, 1, 'compi', [1, 0], ...
                'type', 'rate', 'val', -0.01*meter^3/second);
W.components = info.initial;

% Schedule and initial state
schedule = simpleSchedule(rampupTimesteps(year, 30*day), 'W', W);
state0 = initCompositionalState(model, info.pressure, info.T, [1, 0], info.initial);

% Setup for packed simulations
nm = 'depletion_spe5';
packer = @(model, name, varargin) ...
    packSimulationProblem(state0, model, schedule, nm, 'name', name, varargin{:});

%% First alternative: single separator
% The separator will flash the surface streams at the surface conditions of
% 300 K and one atmosphere pressure, and use the result to determine the
% appropriate volume to extract at reservoir conditions to meet the target
% depletion rate.
model_sep = model;          
s  = EOSSeparator('pressure', 1*atm, 'T', 300); % Set conditions for surface
sg = SeparatorGroup(s);                         % Group = single separator
sg.mode = 'moles';                              % Use mole mode 
model_sep.FacilityModel.SeparatorGroup = sg;    % Connect to reservoir model

%% Second alternative: a tree of separators
% In the second case, the produced reservoir fluids are assumed to undergo
% a more complex separation process before reaching the final stock-tank
% conditions at the surface. The tree of separators is modelled as a
% directed graph using the "dest" array, in which the first column contains
% the liquid targets and the second column the vapor targets.  Each
% separator has a liquid and a vapor target. A target is either a positive
% index indicating the next step of separation or a zero to signify that
% the output stream should be sent to the surface separator (tank) that
% flashes the final gas and liquid streams at surface conditions.
model_multisep = model;

p = [200, 175, 50, 10]*barsa; % p for each sep
T = [423, 400, 350, 300];     % T for each sep
dest =  [2, 3; ... % Send liquid to 2, vapor to 3
         0, 4; ... % Send liquid to final, vapor to 4
         4, 0; ... % Send liquid to separator 4, vapor to tank
         0, 0];    % Send both to tank
sg = SeparatorGroup(s, p, T, dest);
sg.mode = 'moles';
model_multisep.FacilityModel.SeparatorGroup = sg;

%% Simulate the three alternative models
problem      = packer(model, 'Overall-NoSep');
problem_sep  = packer(model_sep, 'Overall-SingleSep');
problem_msep = packer(model_multisep, 'Overall-MultiSep');
problems     = {problem, problem_msep, problem_sep};

simulatePackedProblem(problems);

%% Plot the results
[ws, states, reports, names, T] = getMultiplePackedSimulatorOutputs(problems);
plotWellSols(ws, T, 'datasetnames', names)

%% Plot surface production
nd = numel(ws);
ns = numel(schedule.step.val);
qos = zeros(ns, nd);
qgs = zeros(ns, nd);
for i = 1:nd
    qos(:, i) = getWellOutput(ws{i}, 'qOs');
    qgs(:, i) = getWellOutput(ws{i}, 'qGs');
end
figure
subplot(2,1,1); 
plot(qos,'LineWidth',2); legend(names); title('Liquid rates');
subplot(2,1,2); 
plot(qgs,'LineWidth',2); legend(names); title('Vapor rates');

%% Plot gas saturation and pressure in well cell
lw = 1.5;
c = W.cells;
t = cumsum(schedule.step.val);
figure; hold on
for i = 1:nd
    p = cellfun(@(x) x.pressure(c), states{i});
    sg = cellfun(@(x) x.s(c, 2), states{i});
    yyaxis right
    plot(t/day, sg, 'linewidth', lw);
    ylabel('Gas saturation');
    yyaxis left
    plot(t/day, p/barsa, 'linewidth', lw);
    ylabel('Pressure [bar]')
end
xlabel('Time [days]'); set(gca,'XLim',[0 year/day]);
legend(names)

%%
mrstModule add mrst-gui
figure;
plotToolbar(G, states{2}); view(3);

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