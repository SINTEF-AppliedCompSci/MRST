%% Gravity segregation using two phase AD solvers
mrstModule add ad-core ad-blackoil ad-props ad-fi mrst-gui
[schedule, model, state0] = getBenchmarkAD('spe1');

% lim = 5;
% schedule.step.val = schedule.step.val(1:lim);
% schedule.step.control = schedule.step.control(1:lim);

%% Simulate the problem
close all
fn = getPlotAfterStep(state0, model, schedule, 'plotWell', true, 'plotReservoir', true);

%%
[wellSols, states, report] = simulateScheduleAD(state0, model, schedule, ...
    'Verbose', true, 'afterStepFn', fn);
%%
close all
