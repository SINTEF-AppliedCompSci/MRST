%% Injection of a highly mobile fluid into a fractured medium
% This example is a simple incompressible, two-phase analogue of the
% examples studied in Section 5.2 of Moyner & Tchelepi, SPE J, 23(6), 2018,
% doi: 10.2118/182679-PA. We inject a fluid of density 300 kg/m^3 and
% viscosity 0.3 cP into another fluid of density 1000 kg/m^3 and viscosity
% 1 cP that fills a layered medium containing thirteen high-permeability
% fractures.
%
% The example assumes that you have METIS installed on your system. This is
% relatively simple on a Linux system; see callMetisMatrix function for
% details how to set it up. On Windows, we do unfortunately not yet have a
% general receipe of how to set up METIS.
mrstModule add coarsegrid incomp msrsb book ad-core

%%  Set up the fine-scale model
gravity reset off
rdim = [1000 500];
pth = fullfile(getDatasetPath('MSFractures'), 'setup_fracture.mat');
load(pth);
G.nodes.coords = bsxfun(@times, G.nodes.coords, rdim);
G = computeGeometry(G);
perm(G.cells.tag)=10000;
rock = makeRock(G, perm*milli*darcy, 0.3);
hT   = computeTrans(G, rock);

figure
plotCellData(G, log10(perm),'EdgeAlpha',.1);
mrstColorbar(perm,'south',true); axis tight

fluid = initSimpleFluid('mu', [.3, 1]*centi*poise, ...
    'n', [2, 2], 'rho', [300, 1000]*kilo*gram/meter.^3);

state0 = initResSol(G, 300*barsa, [0, 1]);

icell = findEnclosingCell(G, [0.025, 0.05].*rdim);
pcell = findEnclosingCell(G, [0.975, 0.95].*rdim);
pv    = sum(poreVolume(G, rock));
time  = 10*year;

W  = addWell([], G, rock, icell, 'type', 'rate', ...
             'val',  pv/time, 'compi', [1 0], 'Name', 'I');
W  = addWell(W, G, rock, pcell, 'type', 'rate', ...
             'val', -pv/time, 'compi', [1 0], 'Name', 'P');
plotGrid(G,[icell pcell],'FaceColor','w');

%% Set up two different multiscale configurations
% Rectangular partition
pr     = sampleFromBox(G, reshape(1:100,[10 10]));
pr     = processPartition(G, compressPartition(pr));
CG1    = generateCoarseGrid(G, compressPartition(pr));
CG1    = coarsenGeometry(CG1);
CG1    = storeInteractionRegion(CG1);
basis1 = getMultiscaleBasis(CG1, getIncomp1PhMatrix(G, hT), 'type', 'msrsb');

% METIS partition adapting to fractures
% 
T   = 1 ./ accumarray(G.cells.faces(:,1), 1 ./ hT, [G.faces.num, 1]);
tag = [0; G.cells.tag];
N   = G.faces.neighbors + 1;
T(tag(N(:, 1)) ~= tag(N(:, 2))) = 0;
pa     = partitionMETIS(G, T, 90);
CG2    = generateCoarseGrid(G, compressPartition(pa));
CG2    = coarsenGeometry(CG2);
CG2    = storeInteractionRegion(CG2);
basis2 = getMultiscaleBasis(CG2, getIncomp1PhMatrix(G, hT), 'type', 'msrsb');

% Alternative adapted approach
% pr = sampleFromBox(G, reshape(1:64,[8 8]));
% pf = zeros(G.cells.num,1); pf(G.cells.tag>0)=1;
% pa = processPartition(G, compressPartition(pf*64+pr));
% pa = mergeBlocks(pa, G, ones(G.cells.num,1), ones(G.cells.num,1), 10);

%% Set up solvers
dt = rampupTimesteps(time, time/100, 5);

psolver   = @(u)     incompTPFA(u, G, hT, fluid, 'W', W);
tsolver   = @(u, dt) implicitTransport(u, G, dt, rock, fluid, 'W', W);
mssolver1 = @(u)     incompMultiscale(u, CG1, hT, fluid, basis1, 'W', W);
mssolver2 = @(u)     incompMultiscale(u, CG2, hT, fluid, basis2, 'W', W);

[rstates, rws, mstates1, mws1, mstates2, mws2] = deal(cell(numel(dt)+1,1));
rstates{1}  = psolver(state0);    rws{1} = getWellSol(W, rstates{1},  fluid);
mstates1{1} = mssolver1(state0); mws1{1} = getWellSol(W, mstates1{1}, fluid);
mstates2{1} = mssolver2(state0); mws2{1} = getWellSol(W, mstates2{1}, fluid);

%% Plot the initial saturation distribution
%figure
subplot(1,3,1)
hp(1) = plotCellData(G, rstates{1}.s(:,1),'EdgeColor','none'); 
view(-90,90), axis tight, caxis([0 1]), colormap(flipud(winter))
plotGrid(G,[icell pcell],'FaceColor','r','EdgeColor','r');
title('Fine scale','FontSize',12,'FontWeight','normal');

subplot(1,3,2)
hp(2) = plotCellData(G, mstates2{1}.s(:,1),'EdgeColor','none'); 
view(-90,90), axis tight, caxis([0 1]), colormap(flipud(winter(10).^1.5))
plotGrid(G,[icell pcell],'FaceColor','r','EdgeColor','r');
plotFaces(CG1,1:CG1.faces.num,'EdgeColor','k','FaceColor','none');
title('MsRSB','FontSize',12,'FontWeight','normal');

subplot(1,3,3)
hp(3) = plotCellData(G, mstates2{1}.s(:,1),'EdgeColor','none'); 
view(-90,90), axis tight, caxis([0 1]), colormap(flipud(winter(10).^1.5))
plotGrid(G,[icell pcell],'FaceColor','r','EdgeColor','r');
plotFaces(CG2,1:CG2.faces.num,'EdgeColor','k','FaceColor','none');
title('MsRSB','FontSize',12,'FontWeight','normal');

%% Run the simulation
for i = 1:numel(dt)
    state         = psolver(rstates{i});
    rstates{i+1}  = tsolver(state, dt(i));
    rws{i+1}      = getWellSol(W, rstates{i+1}, fluid);
    hp(1).CData   = rstates{i+1}.s(:,1);
    
    state         = mssolver1(mstates1{i});
    mstates1{i+1} = tsolver(state, dt(i));
    mws1{i+1}     = getWellSol(W, mstates1{i+1}, fluid);
    hp(2).CData   = mstates1{i+1}.s(:,1);

    state         = mssolver2(mstates2{i});
    mstates2{i+1} = tsolver(state, dt(i));
    mws2{i+1}     = getWellSol(W, mstates2{i+1}, fluid);
    hp(3).CData   = mstates2{i+1}.s(:,1);

    drawnow
end
plotWellSols({rws, mws1, mws2}, [0; cumsum(dt)], 'field', 'wcut', ...
    'datasetnames', {'fine scale','MsRSB rectangular','MsRSB adapted'});

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
