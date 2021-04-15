%% Computation of Adjoints for Polymer and Waterflooding
% In this example, we demonstrate how one can use an adjoint simulation to
% compute gradients (sensitivities). To illustrate this, we compare water
% and polymer flooding in a simple box-shaped reservoir represented on a
% uniform 31×31×3 Cartesian grid with homogeneous rock properties. The well
% pattern consists of a producer in the center of the reservoir and two
% injectors in the northeast and southwest corners; all wells are completed
% in the top layer. Somewhat unrealistic, the well schedule consists of a
% period of injection with polymer, followed by a water flooding phase
% without polymer. Finally, the water rate is reduced for the final time
% steps. For comparison, we also simulate a case with pure waterflooding.
%
% The implementation and computational strategy used to obtain adjoints is
% independent of the actual system considered and a similar approach can be
% used for other systems and setups with minor (and obvious) modifications.
%

mrstModule add ad-core ad-blackoil ad-eor ad-props deckformat optimization

%% Read input file and setup basic structures
% The computational setup is described in an input file on the ECLIPSE
% input format. We start by reading this file and extracting the necessary
% structures to represent grid, petrophysics, and fluid properties.

current_dir = fileparts(mfilename('fullpath'));
if mrstIsLiveEditorDir(current_dir)
   % Running in "Live Editor" cell mode.  Fall back to expected location.

   current_dir = fullfile(mrstPath('ad-eor'), 'examples', 'polymer');
end
fn  = fullfile(current_dir, 'POLYMER.DATA');

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

%% Set up initial reservoir state
% Initially, the reservoir is filled by oil on top and water at the bottom.
% The horizontal oil-water contact passes through the second grid layer and
% hence we get different saturations in each of the three grid layers.

ijk = gridLogicalIndices(G);
state0 = initResSol(G, 300*barsa, [ .9, .1]);
state0.s(ijk{3} == 1, 2) = .9;
state0.s(ijk{3} == 2, 2) = .8;

% Enforce s_w + s_o = 1;
state0.s(:,1) = 1 - state0.s(:,2);

% Add zero polymer concentration to the state.
state0.cp    = zeros(G.cells.num, 1);
state0.cpmax = zeros(G.cells.num, 1);

clf, title('Oil concentration')
plotCellData(G, state0.s(:,2));
plotGrid(G, 'facec', 'none')
axis tight off, view(70, 30), colorbar;

%% Plot polymer properties
% In the basic setup, water is much more mobile than the oil we are trying
% to displace, and as a result, a pure waterflooding will not be very
% efficient. To improve the oil recovery, we will add polymer to the
% injected water phase. This will increase the viscosity of the water
% phase and hence make the water much less mobile, which in turn will
% increase the local displacement efficiency and also the volumetric sweep.
% To see the effect of polymer on the water viscosity, we plot the visocity
% multiplier, which is based on tabulated data given in the input file
dc = 0 : .1 : fluid.cpmax;
plot(dc, fluid.muWMult(dc),'LineWidth',2)
title('muW Multiplier')
xlabel('Polymer concentration [kg/m^3]')


%% Set up systems
% To quantify the effect of adding the polymer to the injected water, we
% will solve the same system both with and without polymer. This is done by
% creating an Oil/Water/Polymer system and an Oil/Water system. Note that
% since the data file already contains polymer as an active phase, we do not
% need to pass initADISystem anything other than the deck.

modelPolymer = OilWaterPolymerModel(G, rock, fluid, 'inputdata', deck);
modelOW = TwoPhaseOilWaterModel(G, rock, fluid, 'inputdata', deck);

% Convert the deck schedule into a MRST schedule by parsing the wells
schedule = convertDeckScheduleToMRST(modelPolymer, deck);


%% Run the schedule
% Once a model has been created it is trivial to run the schedule. Any
% options such as maximum non-linear iterations and tolerance can be set in
% the model and the NonLinearSolver class (optional argument to
% simulateScheduleAD).

[wellSolsPolymer, statesPolymer] = simulateScheduleAD(state0, modelPolymer, ...
                                                  schedule);
%%
[wellSolsOW, statesOW] = simulateScheduleAD(state0, modelOW, schedule);

plotWellSols({wellSolsPolymer; wellSolsOW},'datasetnames',{'Polymer','Water'});

%% Objective functions
% Create objective functions for the different systems. We set up
% approximate prices in USD for both the oil and the injection cost of the
% different phases. The polymer injection cost is per kg injected.
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

objectiveInexpensivePolymer = objectivePolymer( 1.0);
objectiveRegularPolymer     = objectivePolymer( 5.0);
objectiveExpensivePolymer   = objectivePolymer(15.0);

%% Plot accumulated present value
% The objective function is now the net present value of the reservoir at
% each time step, i.e. the cost of doing that timestep. However, the most
% interesting value is here the accumulated net present value, as it will
% show us the profit over the lifetime of the reservoir. We plot the three
% different polymer costs, as well as the total profit without polymer
% injection.
%
% While polymer is injected, the net present value may be lower than
% without polymer as there is an increased cost. Once the polymer injection
% phase is over, however, we hope to reap the benefits and get an increased
% oil output resulting in a bigger total value for the reservoir lifetime.

figure(gcf); clf;
cumt = cumsum(schedule.step.val);
v = @(value) cumsum([value{:}]);

plot(convertTo(cumt, year), ...
     convertTo([v(objectiveOW)            ; ...
                v(objectiveInexpensivePolymer)  ; ...
                v(objectiveRegularPolymer); ...
                v(objectiveExpensivePolymer)], 1e6),'LineWidth',2)

legend('Without polymer'       , ...
       'With inexpensive polymer'    , ...
       'With regular polymer'  , ...
       'With expensive polymer', 'Location', 'NorthWest')
title('Net present value')
ylabel('Million USD')
xlabel('Years')
axis tight

%% Compute gradients using the adjoint formulation
% We pass a function handle to the polymer equations and calculate the
% gradient with regards to our control variables. The control variables are
% defined as the last two variables, i.e., well closure (rate/BHP) and
% polymer injection rate.

adjointGradient = ...
    computeGradientAdjointAD(state0, statesPolymer, modelPolymer, ...
                             schedule, objectivePolymerAdjoint);


%% Plot the gradients
% Plot the polymer and well gradients. Note that the second injection well,
% which injects less polymer, should be matched with the first to maximize
% the value. If we were to employ this example in an optimization loop
% using some external algorithm, we could optimize the polymer injection
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

%% Animation of simulation results
% We animate the computed solutions. The animation shows all cells with
% water saturation exceeding the residual saturation, as well as a pie
% chart that shows the oil/water ratio in the producer with and
% without polymer injection in a single chart. 
W = schedule.control(1).W;
ph = nan;
figure(gcf); clf
nDigits = floor(log10(numel(statesPolymer) - 1)) + 1;
subplot(1,3,[1 2])
plotGrid(G, 'facea', 0,'edgea', .05, 'edgec', 'k'); plotWell(G,W);
axis tight off, view(6,60), hs = [];
for i = 1 : numel(statesPolymer) - 1
    injp  = wellSolsPolymer{i}(3);
    injow = wellSolsOW{i}(3);
    state = statesPolymer{i + 1};

    subplot(1, 3, 3), cla
    rates = injp.sign .* [injow.qWs, injow.qOs injp.qWs, injp.qOs];
    
    if ishandle(ph)
        delete(ph)
        ph = pie(rates./sum(rates));
    else
        ph = pie(rates./sum(rates));
        legend('Water (No polymer)', 'Oil (No polymer)', ...
               'Water (With polymer', 'Oil (Polymer)',   ...
               'Location', 'SouthOutside')
        title('Producer OW ratio'),
    end
    subplot(1, 3, [1 2]), delete(hs);
    hs = plotCellData(G,state.s(:,1),state.s(:,1)>.101);
    title(sprintf('Step %0*d (%s)', nDigits, i, formatTimeRange(cumt(i))));
    drawnow;
end

%%
% If you have a powerful computer, you can replace the cell plot by a
% simple volume renderer by uncommenting the corresponding lines.
%{
clf,
plotGrid(G,'facea', 0,'edgea', .05, 'edgec', 'k'); plotWell(G,W);
axis tight off, view(6,60)
for i=1:numel(statesPolymer)-1
    cla
    state = statesPolymer{i + 1};
    plotGrid(G, 'facea', 0,'edgea', .05, 'edgec', 'k');
    plotGridVolumes(G, state.s(:,2), 'cmap', @copper, 'N', 10)
    plotGridVolumes(G, state.s(:,1), 'cmap', @winter, 'N', 10)
    plotGridVolumes(G, state.c,      'cmap', @autumn, 'N', 10)
    plotWell(G,W);
    title(sprintf('Step %0*d (%s)', nDigits, i, formatTimeRange(cumt(i))));
    drawnow
end
%}
    
%% Plot the accumulated water and oil production for both cases
% We concat the well solutions and plot the accumulated producer rates for
% both the polymer and the non-polymer run. The result shows that the
% polymer injection both gives more produced oil and less produced water
% and hence is a feasible stategy in this particular example.
figure(gcf); clf;
wspoly = vertcat(wellSolsPolymer{:});
wsow = vertcat(wellSolsOW{:});

data = -([ [wsow(:,3).qWs]   ; ...
           [wspoly(:,3).qWs] ; ...
           [wsow(:,3).qOs]   ; ...
           [wspoly(:,3).qOs] ] .');
data = bsxfun(@times, data, schedule.step.val);

plot(convertTo(cumt, year), cumsum(convertTo(data, stb)),'LineWidth',2);
legend({'Water without polymer', 'Water with polymer', ...
        'Oil without polymer', 'Oil with polymer'}, 'Location', 'NorthEast');
ylabel('stb');
xlabel('Years');

%% Copyright notice
%
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
