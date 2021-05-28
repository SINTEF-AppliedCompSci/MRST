%% Effect of polymer retention/perm reduction on macroscopic displacement
% The script runs three simulations on a vertical cross-section from the
% SPE10 benchmark to illustrate the effect that polymer retention and
% reduced permeability has on flow conformance. The simulations are:
%
% - Pure water injection
% - Polymer flooding without adsorption
% - Polymer flooding with adsorption and strong permeability reduction
clc
mrstModule add ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui spe10

%% Load data from SPE10
gravity reset on;
rock  = getSPE10rock(1,1:200,1:20);
rock.poro(rock.poro==0) = 1e-4;
[I,J] = deal(repmat([1 200],20,1),1:20);
inx   = sub2ind([200 20],I(:),[J'; J']);
rock.poro(inx)   = 0.2;
rock.perm(inx,1) = 100*milli*darcy;

%% Run case with polymer retention
bookdir = getDatasetPath('eor_book_ii');
fn      = fullfile(bookdir,'macrodisp','vertSPE10','crossect.DATA');
deck    = readEclipseDeck(fn);
deck    = convertDeckUnits(deck);
deck.GRID.PORO  = rock.poro;
deck.GRID.PERMX = rock.perm(:,1);
deck.GRID.PERMY = rock.perm(:,2);
deck.GRID.PERMZ = rock.perm(:,3);
[state0, model, schedule] = initEclipseProblemAD(deck);
model.AutoDiffBackend = SparseAutoDiffBackend;
model.usingShear      = false;
[ws, states, rep]     = simulateScheduleAD(state0, model, schedule);

%% Run case without polymer retention
bookdir = getDatasetPath('eor_book_ii');
fn      = fullfile(bookdir,'macrodisp','vertSPE10','crossect_noads.DATA');
deck    = readEclipseDeck(fn);
deck    = convertDeckUnits(deck);
deck.GRID.PORO  = rock.poro;
deck.GRID.PERMX = rock.perm(:,1);
deck.GRID.PERMY = rock.perm(:,2);
deck.GRID.PERMZ = rock.perm(:,3);
[state0, modelNA, scheduleNA] = initEclipseProblemAD(deck);
modelNA.AutoDiffBackend = SparseAutoDiffBackend;
modelNA.usingShear      = false;
[wsNA, statesNA, repNA] = simulateScheduleAD(state0, modelNA, scheduleNA);


%% Run waterflooding
schedule.control.W(1).cp = 0;
[wsW, statesW, repW]  = simulateScheduleAD(state0, model, schedule);

%%
figure,
nstep = 50;
subplot(2,2,1)
plotCellData(model.G, states{nstep}.s(:,1),'EdgeColor','none')
view(0,0), axis tight, colormap(flipud(winter)); axis off
subplot(2,2,2)
plotCellData(model.G, statesNA{nstep}.s(:,1),'EdgeColor','none')
view(0,0), axis tight, colormap(flipud(winter)); axis off
subplot(2,2,4)
plotCellData(model.G, statesW{nstep}.s(:,1),'EdgeColor','none')
view(0,0), axis tight, colormap(flipud(winter)); axis off

bx = subplot(2,2,3); cla
model = model.validateModel();
kx    = reshape(model.rock.perm(:,1),model.G.cartDims([1 3]))';
rrf   = model.getProp(states{nstep},'PolymerPermReduction');
rrf   = reshape(rrf, model.G.cartDims([1 3]))';
image(log10(kx./rrf), 'CDataMapping','scaled');
colormap(bx,.85*flipud(pink)+.15);
cax=caxis();
hold on, contour(1:200,1:20,rrf,1,'LineWidth',1), hold off
caxis(cax); axis off

axlim = axis(bx);
ax    = axes('Position',bx.Position);
ha    = image(ax, rrf, 'CDataMapping','scaled');
set(ax,'YDir','reverse')
ha.AlphaData = 0.2*(rrf>1.01);
colormap(ax,rot90(autumn,2));
axis(ax,axlim); axis off

subplot(2,2,1);
set(gca,'Position',get(gca,'Position')-[0.00 .1 0 0]);
annotation(gcf,'textbox',[0.36 0.723 0.1 0.097],'Color',[.3 .3 .3],...
    'String',{'Polymer flood','with retention'}, 'LineStyle','none',...
    'FitBoxToText','off');

subplot(2,2,2);
set(gca,'Position',get(gca,'Position')-[0.08 .1 0 0]);
annotation(gcf,'textbox',[0.719 0.723 0.11 0.097],'Color',[.3 .3 .3],...
    'String',{'Polymer flood','without retention'}, 'LineStyle','none',...
    'FitBoxToText','off');

subplot(2,2,4); set(gca,'Position',get(gca,'Position')-[0.08 .0 0 0]);
annotation(gcf,'textbox',[0.719 0.35 0.1 0.097],'Color',[.3 .3 .3],...
    'String',{'Waterflood'}, 'LineStyle','none',...
    'FitBoxToText','off');

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
