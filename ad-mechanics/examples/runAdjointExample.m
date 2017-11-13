%% Computation of Adjoints for Lift values
% In this example, we demonstrate how one can use an adjoint simulation to
% compute gradients (sensitivities). 
%
% The implementation and computational strategy used to obtain adjoints is
% independent of the actual system considered and a similar approach can be
% used for other systems and setups with minor (and obvious) modifications.
%

mrstModule add ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui

%% Setup case


% blackoil

% 2D cartesian domain, try some irregular geometry from run2Dcase

% Same fluid as 2D case (black oil)

%% Setup material parameters for Biot and mechanics

% same as 2D case
E          = 1 * giga * Pascal; % Young's module
nu         = 0.3;               % Poisson's ratio
alpha      = 1;                 % Biot's coefficient
Ev         = repmat(E, G.cells.num, 1);
nuv        = repmat(nu, G.cells.num, 1);
rock.alpha = repmat(alpha, G.cells.num, 1);

%% Setup boundary conditions for mechanics (no displacement)

% fixed bottom, left and right. Fixed at the top.
% todo ...

%% Set up initial reservoir state
clear initState;
initState.pressure = pRef*ones(G.cells.num, 1);
switch opt.fluid_model
  case 'blackoil'
    init_sat = [0, 1, 0];
    initState.rs  = 0.5*fluid.rsSat(initState.pressure);
  case 'oil water'
    init_sat = [0, 1];
  case 'water'
    init_sat = [1];
  otherwise
    error('fluid_model not recognized.')
end
initState.s = ones(G.cells.num, 1)*init_sat;
initState.xd = zeros(nnz(~model.mechModel.operators.isdirdofs), 1);
initState = addDerivedQuantities(model.mechModel, initState);




%% Setup load for mechanics

% In this example we do not impose any volumetric force
loadfun = @(x) (0*x);

%% Gather all the mechanical parameters in a struct

mech = struct('Ev', Ev, 'nuv', nuv, 'el_bc', el_bc, 'load', loadfun);


%% Setup model

model = MechBlackOilModel(G, rock, fluid, mech, 'verbose', true);

%% Setup wells

% 2 injection wells on the side, one production well in middle.

%% Setup a schedule

% constant injection rate.
% todo ...


%% Run the schedule

[wellSols, states] = simulateScheduleAD(state0, model, schedule);


%% Objective functions

% We set up the objective function as beeing the final uplift in one point above the
% production well.

% todo ...

% Create objective functions for the different systems. We set up
% approximate prices in USD for both the oil and the injection cost of the
% different phases. The polymer injection cost is per kg injected.
prices = {'OilPrice',            100  , ...
          'WaterProductionCost',   1  , ...
          'WaterInjectionCost',    0.1, ...
          'DiscountFactor',        0.1 };

objectivePolymerAdjoint = @(tstep) NPVOWPolymer(G, wellSolsPolymer, schedule, ...
                                                'ComputePartials', true, ...
                                                'tStep', tstep, prices{:});

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
for i = 1 : numel(statesPolymer) - 1,
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

% <html>
% <p><font size="-1">
% Copyright 2009-2017 SINTEF ICT, Applied Mathematics.
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
