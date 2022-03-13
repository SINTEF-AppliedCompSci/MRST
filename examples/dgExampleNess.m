%% Simulation of an inverted five-spot pattern with dG
% In this example, we simulate water injection in an inverted five-spot
% pattern posed on a subset of SPE10 using dG(0) and dG(1) visualize the
% results

%% Add modules
mrstModule add dg spe10 ad-props ad-core ad-blackoil sequential

%% Set up fine-scale model
% We extract part of layer 51 of SPE10 2
[state0, imodel, schedule] = ...
    setupSPE10_AD('layers', 51, 'J',1:110,'make2D', true, 'T', 3*year, 'dt', 20*day);
G  = computeCellDimensions(imodel.G);
schedule.control.W(3).cells = 6550;
schedule.control.W(4).cells = 6578;
schedule.control.W(5).cells = 3811;

%% Set base model
model   = GenericBlackOilModel(G, imodel.rock, imodel.fluid, 'gas', false);
pmodel  = PressureModel(model); % Pressure model
tmodel0 = TransportModelDG(model, 'degree', 0);
model0  = SequentialPressureTransportModel(pmodel, tmodel0,'parentModel', model);
model0.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true);
[ws0, states0, reports0] = simulateScheduleAD(state0, model0, schedule);

%% Simulate with dG(1)
tmodel1 = TransportModelDG(model, 'degree', 1, 'dsMaxAbsDiv', 2);
model1  = SequentialPressureTransportModel(pmodel, tmodel1,'parentModel', model);
model1.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true);
model1.transportNonLinearSolver = NonLinearSolver('useLinesearch', true);
[ws1, states1, reports1] = simulateScheduleAD(state0, model1, schedule);

%% Make animation of the test case
% To produce nice plots, we use three axes on top of each other. The first
% plots the permeability using a pink colormap, the second plots cells with
% water saturation exceeding 0.2 using a green-blue colormap, whereas the
% third only shows the well positions.
figure('Position', [0,0,900,400])
xw = model.G.cells.centroids(vertcat(schedule.control.W.cells),1:2);
bx1 = subplot(1,2,1); 
plotCellData(G,log10(model.rock.perm(:,1)),'EdgeColor','none'); axis tight off
colormap(bx1,flipud(pink));
axlim = axis(bx1);
ax1 = axes('Position',bx1.Position);
colormap(ax1,flipud(winter));
cx1 = axes('Position',bx1.Position);
plot(xw(:,1),xw(:,2),'.r','MarkerSize',18);
axis(cx1,axlim); axis off

bx2 = subplot(1,2,2); 
plotCellData(G,log10(model.rock.perm(:,1)),'EdgeColor','none'); axis tight off
colormap(bx2,flipud(pink));
ax2 = axes('Position',bx2.Position);
colormap(ax2,flipud(winter));
cx2 = axes('Position',bx2.Position);
plot(xw(:,1),xw(:,2),'.r','MarkerSize',18);
axis(cx2,axlim); axis off


[h1,h2]=deal([]);
for i=1:numel(states0)
    delete([h1 h2]);
    set(gcf,'CurrentAxes',ax1);
    h1 = plotCellData(G,states0{i}.s(:,1),states0{i}.s(:,1)>.205,'EdgeColor','none');
    axis(ax1,axlim), axis off
    set(gcf,'CurrentAxes',ax2);
    h2 = plotCellData(G,states1{i}.s(:,1),states1{i}.s(:,1)>.205,'EdgeColor','none');
    axis(ax2,axlim), axis off
    drawnow;
    pause(.2)
end

%% Plot production curves
plotWellSols({ws0, ws1},reports0.ReservoirTime,...
    'datasetnames',{'dG(0)','dG(1)'}, 'field', 'qWs');

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
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
