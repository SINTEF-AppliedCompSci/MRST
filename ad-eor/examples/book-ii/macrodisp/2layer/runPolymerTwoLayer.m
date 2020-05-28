%% Conformance in a layered system
% This script demonstrates how injecting polymer gives improved flow
% conformance in a layered system by improving the mobility ratio between
% displacing and displaced fluids.
clc
mrstModule add ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui

%% Run polymer flooding
gravity reset on;
bookdir = getDatasetPath('eor_book_ii');
fn      = fullfile(bookdir,'macrodisp','2layer','polymerTwoLayer.DATA');
[state0, model, schedule, nonlinear] = initEclipseProblemAD(fn);
model.usingShear = false;

fn = getPlotAfterStep(state0, model, schedule, 'plotWell', true, ...
    'plotReservoir', true, 'field', 's:1', 'view', [0, 0]);
[wsP, statesP, reportP] = ...
    simulateScheduleAD(state0, model, schedule, ...
    'NonLinearSolver', nonlinear, 'afterStepFn', fn);

%% Run waterflooding
schedule.control.W(1).cp = 0;
fn = getPlotAfterStep(state0, model, schedule, 'plotWell', true, ...
    'plotReservoir', true, 'field', 's:1', 'view', [0, 0]);
[wsW, statesW, reportW] = ...
    simulateScheduleAD(state0, model, schedule, ...
    'NonLinearSolver', nonlinear, 'afterStepFn', fn);

%% Plot solutions
figure('Position',[20 460 1320 280]);
subplot(1,2,1)
plotCellData(model.G, statesW{end}.s(:,1),'EdgeColor','none'); view(0,0)
hold on, plot3([0 60],-[.1 .1],[1005 1005],'--','Color',[.3 .3 .3]); hold off
text(45,-.1,1009,'waterflood','Color',[.3 .3 .3]);
text(6,-.1,1009,'low perm','Color',[.3 .3 .3]);
text(6,-.1,1001,'high perm','Color',[.3 .3 .3]);

subplot(1,2,2)
plotCellData(model.G, statesP{end}.s(:,1),'EdgeColor','none'); view(0,0)
hold on, plot3([0 60],-[.1 .1],[1005 1005],'--','Color',[.3 .3 .3]); hold off
text(45,-.1,1009,'polymer flood','Color',[.3 .3 .3]);

colormap(flipud(winter));
cb=colorbar('Location','EastOutside');
set(cb,'Position',[.925 .2 .02 0.6]);

