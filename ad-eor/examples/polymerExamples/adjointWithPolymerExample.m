%% Read case from file
% This example contains a simple $31\times31\times3$ fine grid containing
% two injectors in opposite corners and one producer in the middle of the
% domain. All wells are completed in the top layers of cells.
%
% The schedule being used contains first a period of injection with
% polymer, followed by a water flooding phase without polymer. Finally, the
% water rate is reduced for the final time steps.
%

try
    require add ad-core ad-blackoil ad-eor ad-props deckformat optimization
catch
    mrstModule add ad-core ad-blackoil ad-eor ad-props deckformat optimization
end

current_dir = fileparts(mfilename('fullpath'));
fn    = fullfile(current_dir, 'POLYMER.DATA');

deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);

G = initEclipseGrid(deck);
G = computeGeometry(G);

rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);

fluid = initDeckADIFluid(deck);

% Oil rel-perm from 2p OW system.
% Needed by equation implementation function 'eqsfiOWExplictWells'.
fluid.krO = fluid.krOW;

gravity on

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

% Add zero polymer concentration to the state.
state0.c    = zeros(G.cells.num, 1);
state0.cmax = zeros(G.cells.num, 1);

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
% $\lambda_w = \frac{K}{mu_w}$, this makes the water much less mobile. As a
% problem with water injection with regards to oil production is that the
% water is much more mobile than the hydrocarbons we are trying to
% displace, injecting a polymer may be beneficial towards oil recovery.
dc = 0:.1:fluid.cmax;
plot(dc, fluid.muWMult(dc))
title('muW Multiplier')
xlabel('Polymer concentration')
ylabel('kg/m^3')


%% Set up systems.
% To quantify the effect of adding the polymer to the injected water, we
% will solve the same system both with and without polymer. This is done by
% creating both a Oil/Water/Polymer system and a Oil/Water system. Note
% that since the data file already contains polymer as an active phase we
% do not need to pass initADISystem anything other than the deck.

modelPolymer = OilWaterPolymerModel(G, rock, fluid, 'inputdata', deck);
modelOW = TwoPhaseOilWaterModel(G, rock, fluid, 'inputdata', deck);

% Convert the deck schedule into a MRST schedule by parsing the wells
schedule = convertDeckScheduleToMRST(G, modelPolymer, rock, deck);


%% Run the schedule
% Once a system has been created it is trivial to run the schedule. Any
% options such as maximum non-linear iterations and tolerance can be set in
% the system struct.

[wellSolsPolymer, statesPolymer] = ...
    simulateScheduleAD(state0, modelPolymer, schedule);

[wellSolsOW, statesOW] = simulateScheduleAD(state0, modelOW, schedule);


%% Objective functions
% Create objective functions for the different systems. We set up
% approximate prices in USD for both the oil price and the injection cost
% of the different phases. The polymer injection cost is per kg injected.
prices = {'OilPrice',            100  , ...
          'WaterProductionCost',   1  , ...
          'WaterInjectionCost',    0.1, ...
          'DiscountFactor',        0.1 };

objectivePolymerAdjoint = ...
    @(tstep) NPVOWPolymer(G, wellSolsPolymer, schedule, ...
                          'ComputePartials', true, 'tStep', tstep, ...
                          prices{:});

% We first calculate the NPV of the pure oil/water solution.
objectiveOW = NPVOW(G, wellSolsOW, schedule, prices{:});

% Calculate the objective function for three different polymer prices
objectivePolymer = ...
    @(polyprice) NPVOWPolymer(G, wellSolsPolymer, schedule, prices{:}, ...
                              'PolymerInjectionCost', polyprice);

objectiveCheapPolymer     = objectivePolymer( 1.0);
objectiveRegularPolymer   = objectivePolymer( 5.0);
objectiveExpensivePolymer = objectivePolymer(15.0);

%% Plot accumulated present value
% In each time step the objective function is now the net present value of
% the reservoir, i.e. the cost of doing that timestep. However, the most
% interesting value is here the accumulated net present value, as it will
% show us the profit for the lifetime of the reservoir. We plot the three
% different polymer cost as well as the total profit without polymer
% injection.
%
% While polymer injection is happening, the polymer value is lower than
% without polymer as there is an increased cost. Once the polymer injection
% phase is over, we reap the benefits and get an increased oil output
% resulting in a bigger total value for the reservoir lifetime.

figure(gcf); clf;

cumt = cumsum(schedule.step.val);

v = @(value) cumsum([value{:}]);


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
% defined as the last two variables, i.e. well closure (rate/BHP) and
% polymer injection rate.

adjointGradient = ...
    computeGradientAdjointAD(state0, statesPolymer, modelPolymer, ...
                             schedule, objectivePolymerAdjoint);


%% Plot the gradients
% Plot the polymer and well gradients. Note that the second injection well
% which injects less polymer should be matched with the first to maximize
% the value. If we were to employ this example in an optimization loop
% using some external algorithm we could optimize the polymer injection
% rate with regards to the total reservoir value.

figure(gcf); clf;
ng = numel(adjointGradient);
for i = 1:ng
    subplot(2,ng,i)
    plot(adjointGradient{i}(1:2), '*');
    title(['Polymer step ' num2str(i) ])
    set(gca, 'xtick', [1;2], 'xlim', [0.5;2.5]);
    axis 'auto y'
    subplot(2,ng,i+3)
    plot(adjointGradient{i}(3:end), '*'); axis tight
    title(['Wells control step ' num2str(i) ])
    set(gca, 'xtick', [1;2;3], 'xlim', [0.5; 3.5]);
    axis 'auto y'
    xlabel('Well #')
end

%% Plot the schedule
% We visualize the schedule and plot both the water, oil and polymer
% concentration using a simple supplied volume renderer. At the same time
% we visualize the sum of the oil/water ratio in the producer for both
% with and without polymer injection in a single pie chart.

W = schedule.control(1).W;
figure(gcf); clf; 
view(10,65);
drawnow;

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

    view(az, el);
    drawnow;
end

%% Plot the accumulated water and oil production for both cases
% We concat the well solutions and plot the accumulated producer rates for
% both the polymer and the non-polymer run. The result shows that

figure(gcf); clf;
wspoly = vertcat(wellSolsPolymer{:});
wsow = vertcat(wellSolsOW{:});

data = -([ [wsow(:,3).qWs]   ; ...
           [wspoly(:,3).qWs] ; ...
           [wsow(:,3).qOs]   ; ...
           [wspoly(:,3).qOs] ] .');
data = bsxfun(@times, data, schedule.step.val);

plot(convertTo(cumt, year), convertTo(data, stb));
legend({'Water without polymer', 'Water with polymer', ...
        'Oil without polymer', 'Oil with polymer'}, 'Location', 'NorthEast');
ylabel('stb');
xlabel('Years');

