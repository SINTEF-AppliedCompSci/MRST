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

%% Copyright Notice
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
