%% NPV - analysis for a simple 2D model  
% In this example we consider a simple 2D oil-water model. The main purpose
% is a comparison between a base simulation and an optimal NPV simulation
% as obtained by the script 'optimizeModel2D'. The optimization will run
% only the first time this script is called, as the optimal schedule will
% be saved and loaded at next execution. The script is intended for
% MATLAB-Publish, and it is advised to set *Max # of output lines* in 
% *Publish settings* to 1.

% Load required modules 
mrstModule add ad-core ad-blackoil ad-props optimization spe10

% Start by running script in utils-directory to set up example model
setupModel2D
% Plot grid with x-permeability and wells:
close all
figure, plotCellData(G, log(rock.perm(:,1)));
plotWell(G, W, 'Fontsize', 16, 'Color', 'k');
axis off, axis equal, axis tight, camproj perspective, view([-1 -2 1.5]);

%% Base and optimal schedules

% Check if we need to run optimization, if not, load optimal schedule:
try
    load('schedule_opt');
catch
    fprintf('Running optimization, this may take a few minutes ...\n')
    optimizeModel2D;
end
close all
% Display base and optimal schedule:
li = [0 400]/day;        % Injector limits  
lp = [100 250]*barsa;    % Producer limits 
figure,
plotSchedules(schedule, 'singlePlot', true, 'boxConst', [li;li;lp;lp] )
figure,
plotSchedules(schedule_opt, 'singlePlot', true, 'boxConst', [li;li;lp;lp] )

%% Run forward simulation for base and optimal schedule

% Create model-object of class TwoPhaseOilWaterModel
model  = TwoPhaseOilWaterModel(G, rock, fluid);
% Set initial state and run simulations:
state0 = initResSol(G, 200*barsa, [0, 1]);
[wellSols, states]         = simulateScheduleAD(state0, model, schedule);
[wellSols_opt, states_opt] = simulateScheduleAD(state0, model, schedule_opt);
%%
% Plot oil and water production for producer wells for both scenarios:
dtime = schedule.step.val/day;
time  = cumsum(dtime);
qOs     = - getWellOutput(wellSols, 'qOs', 3:4);
qOs_opt = - getWellOutput(wellSols_opt, 'qOs', 3:4);
figure, plot(time, qOs*day, '--', 'LineWidth', 2)
hold on, 
ax = gca; 
if ~isnumeric(ax)
    ax.ColorOrderIndex = 1;
end

plot(time, qOs_opt*day, '-'  ,'LineWidth', 2)
axis([0 600 0 550]), set(gca, 'FontSize', 14)
legend([W(3).name, '-base'] , [W(4).name, '-base'], [W(3).name, '-opt'] , [W(4).name, '-opt'] )
title('Oil production rates')

qWs     = - getWellOutput(wellSols, 'qWs', 3:4);
qWs_opt = - getWellOutput(wellSols_opt, 'qWs', 3:4);
figure, plot(time, qWs*day, '--','LineWidth', 2)
hold on, 
ax = gca; 
if ~isnumeric(ax)
    ax.ColorOrderIndex = 1;
end

plot(time, qWs_opt*day, '-'  ,'LineWidth', 2)
axis([0 600 0 550]), set(gca, 'FontSize', 14)
legend([W(3).name, '-base'] , [W(4).name, '-base'], [W(3).name, '-opt'] , ...
       [W(4).name,  '-opt'], 'Location', 'northwest');
title('Water production rates')

%% Problem economics:
% We here define economic parameters, and analyse net cash flow and NPV.
% We define the _critical water cut_ to be the highest water cut resulting
% in a non-negative net cash-flow under assumption of reservoir volume 
% balance, i.e the water cut for wich we have 
%
% $$r_oq_o^s - r_{wp}q_{wp}^s - r_{wi}q_{wi}^s = 0.$$
% 
% Combining with reservoir volume balance, we obtain:
%
% $$\frac{q_{wp}}{q_o+q_{wp}}\leq\frac{b_or_o-b_wr_{wi}}{b_o(r_o+r_{wp}+r_{wi})-b_wr_{wi}}.$$

% Revenue, costs and discount rate
d   = 0.05;    % yearly discount factor
ro  = 60;      % oil revenue/price ($/stb)
rwp =  6;      % water production handling costs ($/stb)
rwi =  6;      % water injection cost ($/stb) 
npvopts = {'OilPrice',             ro , ...
           'WaterProductionCost', rwp , ...
           'WaterInjectionCost',  rwi , ...
           'DiscountFactor',        d};

% Compute critical water-cut value at 200 bar: 
[bw, bo] = deal(fluid.bW(200*barsa), fluid.bO(200*barsa));  
wcut_crit = (bo*ro-bw*rwi)/(bo*(ro+rwp+rwi)-bw*rwi);
fprintf('Maximal economic water cut: %4.3f\n', wcut_crit)

% Plot water-cut for each producer and critical water-cut:
figure, hold on, plot(time, qWs./(qOs+qWs), 'LineWidth', 2)
axis([0 600 0 1]), set(gca, 'FontSize', 14)
legend(W(3).name, W(4).name, 'Location', 'southeast')
plot([0 600], wcut_crit*[1 1], '--k', 'LineWidth', 2);
title('Water cut - Base')
figure, hold on, plot(time, qWs_opt./(qOs_opt+qWs_opt), 'LineWidth', 2)
axis([0 600 0 1]), set(gca, 'FontSize', 14)
legend(W(3).name, W(4).name, 'Location', 'southeast')
plot([0 600], wcut_crit*[1 1], '--k', 'LineWidth', 2);
title('Water cut - Optimal')

%% Fractional flow contours
% We can also plot fractional flow in the reservoir together with a contour of
% of the critical water cut. We chose the reservoir states at the end of
% each of the four control-steps, and compare _economic_ regions of the
% reservoir:
[~,n] = rlencode(schedule.step.control);
for k = 1:4
    figure,
    fracFlowContours(G, W, states(sum(n(1:k))), fluid, wcut_crit, 'LineWidth', 3, 'color', 'k');
    title(['Base - control step ', num2str(k)])
    figure,
    fracFlowContours(G, W, states_opt(sum(n(1:k))), fluid, wcut_crit, 'LineWidth', 3, 'color', 'k');
    title(['Optimal - control step ', num2str(k)])
end

%% Compute and plot net cashflow and NPV
% Net present value (NPV) is the sum (integral) of discounted net cash
% flows. Hence, for positive cash-flows, NPV is increasing.
close all
vals     = cell2mat(NPVOW(model, states, schedule, npvopts{:}));
vals_opt = cell2mat(NPVOW(model, states_opt, schedule_opt, npvopts{:}));
% Plot discounted net cashflow $/day: 
figure,  plot(time, vals./dtime, '--b','LineWidth', 2);
hold on, plot(time, vals_opt./dtime, '-b','LineWidth', 2);
line([0 600], [0 0], 'color', 'r'), set(gca, 'FontSize', 14)
title('Net cash-flow [$]'), legend('Base', 'Optimal')
% Find index of first occuring time > 10 days, where net cashflow becomes
% negative:
inx = find(and(vals<0, time>10), 1, 'first');

% Plot evolution of NPV and indicate peak value:
npv = cumsum(vals);
figure,  plot(time, cumsum(vals), '--b', 'LineWidth', 2);
hold on, plot(time, cumsum(vals_opt), '-b', 'LineWidth', 2);
plot([1 1]*time(inx), [0 npv(inx)], '--k', 'LineWidth', 2)
set(gca, 'FontSize', 14), title('Evolution of NPV [$]'),
legend('Base', 'Optimal', 'Location', 'northwest')


%% Visualize adjoint gradient for base simulation.
% We compute the gradient by performing an adjoint simulation w.r.t to the 
% NPV-objective. 
objh = @(tstep,model, state) NPVOW(model, states, schedule, 'ComputePartials', true, 'tStep', tstep, npvopts{:});
g     = computeGradientAdjointAD(state0, states, model, schedule, objh);

% We visualize the obtained gradient in a shedule-plot. Actual gradient
% values are scaled for plotting purposes:
figure
plotSchedules(schedule, 'grad', g, 'singlePlot', true, 'boxConst', [li;li;lp;lp] )


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
