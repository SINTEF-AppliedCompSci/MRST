%% 2 phases system with polymer
% In this example, we look at a 2 phases system with polymer. The polymer
% model uses the <http://dx.doi.org/10.2118/3484-PA Todd-Longstaff mixing
% model> and also accounts for dead pore space, adsorption and reduced
% permeability. All these properties are given in the file |poly.inc| which
% follows Eclipse format. Polymer are used to increase the water viscosity
% and establish a more favorable ratio between the oil and water
% mobilities. This example contains a simple $31\times31\times3$ fine grid
% with two injectors in opposite corners and one producer in the middle of
% the domain. All wells are completed in the top layers of cells. We build
% up a polymer plug by first injecting water with polymer (concentration
% $c=c_{max}$ and $c=\frac14c_{max}$ in the first and second injection
% wells, respectively) and, then, water without polymer to maintain the
% pressure gradient in the reservoir.


%% Read case from file
%

% Required modules
mrstModule add deckformat ad-core ad-fi optimization ad-props

current_dir = fileparts(mfilename('fullpath'));
fn    = fullfile(current_dir, 'polymer.data');

deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);

G = initEclipseGrid(deck);
G = computeGeometry(G);

rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);

fluid = initDeckADIFluid(deck);

gravity on


%% Plot wells and permeability
% 

figure(1)
clf;
W = processWells(G, rock, deck.SCHEDULE.control(1));
plotCellData(G, convertTo(rock.perm(:,1), milli*darcy), 'FaceAlpha', .5, ...
            'EdgeAlpha', .3, 'EdgeColor', 'k');
plotWell(G, W);
title('Permeability (mD)')
axis tight;
view(35, 40);
colorbar('SouthOutside');


%% Set up simulation parameters
% We want a layer of oil on top of the reservoir and water on the bottom.
% To do this, we alter the initial state based on the logical height of
% each cell. The resulting oil concentration is then plotted.

ijk = gridLogicalIndices(G);

state0 = initResSol(G, 300*barsa, [ .9, .1]);
state0.s(ijk{3} == 1, 2) = .9;
state0.s(ijk{3} == 2, 2) = .8;

% Enforce s_w + s_o = 1;
state0.s(:,1) = 1 - state0.s(:,2);

% Add zero polymer concentration to the initial state.
state0.c    = zeros(G.cells.num, 1);
state0.cmax = zeros(G.cells.num, 1);

figure(2)
clf
plotCellData(G, state0.s(:,2));
plotGrid(G, 'facec', 'none')
title('Oil concentration')
axis tight off
view(70, 30);
colorbar;

%% Plot polymer properties
% When polymer is added to the water phase, the viscosity of the water
% phase containing polymer is increased. Because mobility is defined as
% $\lambda_w = \frac{K}{mu_w}$, this makes the water less mobile. Polymer
% is used in cases where the water is much more mobile than the
% hydrocarbons we are trying to displace. The viscosity of the polymer
% solution depends on the polymer concentration. In the fully mixed case
% (mixing parameter $\omega=1$), then the effective water viscosity is
% given by $\mu_{w,eff}(c)=\mu_w\mu_c(c)$ where $\mu_c$ is the viscosity
% multiplier and $\mu_w$ the water viscosity without polymer.
figure(3)
clf
dc = 0:.1:fluid.cmax;
plot(dc, fluid.muWMult(dc))
title('viscosity multiplier $\mu_c$')
xlabel('Polymer concentration')
ylabel('kg/m^3')

%% Set up systems.
% To quantify the effect of adding the polymer to the injected water, we
% will solve the same system both with and without polymer. This is done by
% creating both a Oil/Water/Polymer system and a Oil/Water system. Note
% that since the data file already contains polymer as an active phase we
% do not need to pass initADISystem anything other than the deck.

schedule = deck.SCHEDULE;
systemPolymer = initADISystem(deck, G, rock, fluid, 'relaxRelTol', .8);
systemOW =      initADISystem({'Oil', 'Water'}, G, rock, fluid);

systemPolymer.activeComponents %# ok, intentional display
systemOW.activeComponents      %# ok

%% Run the schedule
% Once a system has been created it is trivial to run the schedule. Any
% options such as maximum non-linear iterations and tolerance can be set in
% the system struct.

[wellSolsPolymer, statesPolymer] = ...
   runScheduleADI(state0, G, rock, systemPolymer, schedule);

[wellSolsOW, statesOW] = ...
   runScheduleADI(state0, G, rock, systemOW, schedule);

%% Objective functions
% Create objective functions for the different systems. We set up
% approximate prices in USD for both the oil price and the injection cost
% of the different phases. The polymer injection cost is 0.1 dollar per kg
% injected.

prices = {'OilPrice',            100  , ...
          'WaterProductionCost',   1  , ...
          'WaterInjectionCost',    0.1, ...
          'DiscountFactor',        0.1 };

% We first calculate the NPV of the pure oil/water solution.
objectiveOW = NPVOW(G, wellSolsOW, schedule, prices{:});

% Calculate the objective function for three different polymer prices
objectivePolymer = ...
   @(polyprice) NPVOWPolymer(G, wellSolsPolymer, schedule, prices{:}, ...
                             'PolymerInjectionCost', polyprice);

objectiveCheapPolymer     = objectivePolymer( 1.0);
objectiveRegularPolymer   = objectivePolymer( 5.0);
objectiveExpensivePolymer = objectivePolymer(15.0);

% We define the objective function for the adjoint computation by setting the option
% 'computePartials' to true.
objectivePolymerAdjoint = ...
   @(tstep) NPVOWPolymer(G, wellSolsPolymer, schedule, ...
                         'ComputePartials', true, 'tStep', tstep, ...
                         prices{:}, 'PolymerInjectionCost', 5);



%% Plot accumulated present value
% In each time step the objective function is now the net present value of
% the reservoir, i.e. the cost of doing that timestep. However, the most
% interesting value is here the accumulated net present value, as it will
% show us the profit for the lifetime of the reservoir. We plot the three
% different polymer cost as well as the total profit without polymer
% injection.
%
% While polymer injection is happening, the value with polymer is lower
% because there is an increased cost due to the polymer price. Once the
% polymer injection phase is over, we reap the benefits and get an
% increased oil output resulting in a bigger total value for the reservoir
% lifetime.

cumt = cumsum(schedule.step.val);

v = @(value) cumsum([value{:}]);

figure(4)
clf
plot(convertTo(cumt, year), ...
     convertTo([v(objectiveOW)            ; ...
                v(objectiveCheapPolymer)  ; ...
                v(objectiveRegularPolymer); ...
                v(objectiveExpensivePolymer)], 1e6))

legend('Without polymer'       , ...
       'With cheap polymer'    , ...
       'With regular polymer'  , ...
       'With expensive polymer')

title('Net present value')
ylabel('Million USD')
xlabel('Years')

%% Compute gradient using the adjoint formulation
% We pass a function handle to the polymer equations and calculate the
% gradient with regards to our control variables. The control variables are
% defined as the last two variables, i.e. polymer injection concentration
% and well closure (in this case, all the wells are controled by rate).

ctrl = [6, 7];
adjointGradient = ...
   runAdjointADI(G, rock, fluid, schedule, objectivePolymerAdjoint, ...
                 systemPolymer, 'ForwardStates', statesPolymer,     ...
                 'ControlVariables', ctrl, 'Verbose', true);

%% Plot the gradients
% Plot the polymer and well gradients. Note that the second injection well
% which injects less polymer should be matched with the first to maximize
% the value. If we were to employ this example in an optimization loop
% using some external algorithm we could optimize the polymer injection
% rate with regards to the total reservoir value.

figure(5)
clf
ng = numel(adjointGradient);
for i = 1: ng
    h = subplot('position', [(i-1)/ng, 0, 1/ng, 0.9]);
    set(h, 'outerposition', [(i-1)/ng, 0, 1/ng, 0.9]);
    plot(adjointGradient{i}(1:2), '*');
    title(['step ' num2str(i) ]);
    set(gca, 'xtick', [1; 2], 'xlim', [0.5; 2.5]);
    xlabel('Injection well number');
    ylabel('dollar/(polymer concentration)');
end
h = subplot('position', [0.5, 0.95, 0.01, 0.01]);
set(h,'visible', 'off')
text(0, 0, ['Gradient of total NPV with respect to polymer concentration injection at ' ...
            '' ']each step'], 'fontsize', 14, 'horizontalalignment', 'center', 'fontweight', ...
     'bold');

figure(6)
clf
for i = 1: ng
    h = subplot('position', [(i-1)/ng, 0, 1/ng, 0.9]);
    set(h, 'outerposition', [(i-1)/ng, 0, 1/ng, 0.9]);
    plot(adjointGradient{i}(3:end)*day, '*'); axis tight
    title(['step ' num2str(i) ])
    set(gca, 'xtick', [1; 2; 3], 'xlim', [0.5; 3.5]);
    axis 'auto y'
    ylabel('dollar/(m^3/day)')
    xlabel('Well number')
end
h = subplot('position', [0.5, 0.95, 0.01, 0.01]);
set(h,'visible', 'off')
text(0, 0, ['Gradient of total NPV with respect to injection and production rates at ' ...
            'each step'], 'fontsize', 14, 'horizontalalignment', 'center', 'fontweight', ...
     'bold');


%% Plot the schedule
% We visualize the schedule and plot both the water, oil and polymer
% concentration using a simple supplied volume renderer. At the same time
% we visualize the sum of the oil/water ratio in the producer for both
% with and without polymer injection in a single pie chart.

W = processWells(G, rock, schedule.control(schedule.step.control(1)));
figure(7)
clf
view(10,65)

[az, el] = deal(6, 60);

nDigits = floor(log10(numel(statesPolymer) - 1)) + 1;

for i = 1 : numel(statesPolymer) - 1,
    injp  = wellSolsPolymer{i}(3);
    injow = wellSolsOW{i}(3);

    state = statesPolymer{i + 1};

    subplot(1, 3, 3)
    rates = injp.sign .* [injow.qWs, injow.qOs injp.qWs, injp.qOs];
    pie(rates./sum(rates))
    legend('Water (No polymer)', 'Oil (No polymer)', ...
           'Water (With polymer', 'Oil (Polymer)',   ...
           'Location', 'SouthOutside')

    title('Producer OW ratio')
    subplot(1, 3, [1 2])
    cla
    plotGrid(G, 'facea', 0,'edgea', .05, 'edgec', 'k');
    plotGridVolumes(G, state.s(:,2), 'cmap', @copper, 'N', 10)
    plotGridVolumes(G, state.s(:,1), 'cmap', @winter, 'N', 10)
    plotGridVolumes(G, state.c,      'cmap', @autumn, 'N', 10)
    plotWell(G, W);
    axis tight off

    title(sprintf('Step %0*d (%s)', nDigits, i, formatTimeRange(cumt(i))));

    view(az, el)
    drawnow
end

%% Plot the accumulated water and oil production for both cases
% We concat the well solutions and plot the production rates for
% both the polymer and the non-polymer run. 

figure(8)
clf;
wspoly = vertcat(wellSolsPolymer{:});
wsow = vertcat(wellSolsOW{:});

data = -([ [wsow(:,3).qWs]   ; ...
           [wspoly(:,3).qWs] ; ...
           [wsow(:,3).qOs]   ; ...
           [wspoly(:,3).qOs] ] .');

% We process the values to do bar plots.
cumtbar = cumt(1 : end - 1);
cumtbar = rldecode(cumtbar, 2*ones(numel(cumtbar), 1));
cumtbar = [0;  cumtbar; cumt(end)];
databar = rldecode(data, 2*ones(size(data, 1), 1));

plot(convertTo(cumtbar, year), convertTo(databar, stb));
legend('Water without polymer', 'Water with polymer', ...
       'Oil without polymer', 'Oil with polymer', 'NorthWest')
title('Production rates')
ylabel('stb')
xlabel('Years')

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
