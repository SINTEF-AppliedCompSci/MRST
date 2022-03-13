%% Multiscale solver applied to high-resolution bed model
% Pinchouts will create unstructured non-neighboring connections and hence
% be one of the principal gridding challenges for real-life reservoir
% models. To exemplify, we consider a highly detailed, core-scale model of
% realistic bedding structures. Such models are used as input to derive
% directional permeability for a given lithofacies and identify net pay
% below the level of petrophysical log resolution. The bedding structure
% consists of six different rock types and is realized on a 30x30x333
% corner-point grid. Almost all cells are affected to some degree by pinch:
% Although the model has approximately 300.000 cells initially, over 2/3
% will be inactive due to significant erosion, giving a fine grid consiting
% of approximately 30x30x100 cells. The volumes of the cells, as well as
% the areas of the (vertical) faces, vary almost four orders of magnitude.
% This makes the model challenging to partition.
%
% The model not only has a difficult geometry, but also has a large number
% of low-permeable shale layers pinched between the other high-permeable
% layers. Impermeable regions like this are known to pose monotonicity
% problems for the MsFV method. Hereine, we demonstrate various partition
% strategies and show that MsRSB may also suffer from severe monotonicity
% violations in special cases.

mrstModule add mrst-gui incomp coarsegrid msrsb
gravity off

%% Load model
% This is one of the standard data sets that come along with MRST and can
% be downloaded automatically from our webpages.
pth    = getDatasetPath('bedmodel2');
grdecl = readGRDECL(fullfile(pth, 'BedModel2.grdecl'));
% pth    = getDatasetPath('bedmodels1');
% grdecl = readGRDECL(fullfile(pth,'testModel1.grdecl'));
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);

G = processGRDECL(grdecl);
G = computeGeometry(G);

rock = grdecl2Rock(grdecl, G.cells.indexMap);

%% Plot the petrophysical data 
% The model contains rock type information in the form of the SATNUM
% keyword. Plot the permeability, the facies, as well as the number of
% cells and pore volume per layer in the grid
Kx       = convertTo(rock.perm(:,1),milli*darcy);
facies  = grdecl.SATNUM(G.cells.indexMap);
[~,~,K] = gridLogicalIndices(G);
ncell   = accumarray(K,1);
pvols   = accumarray(K,poreVolume(G,rock));
vols    = accumarray(K,G.cells.volumes);
figure('Position',[300 340 780 420]);
subplot(2,2,1,'FontSize',12)
plotCellData(G,log10(Kx),'EdgeAlpha',.1);
mrstColorbar(Kx,'west',true)
view(30,30); axis tight off;
title('Horizontal permeability','FontSize',10,'FontWeight','normal');

subplot(2,2,3,'FontSize',12)
plotCellData(G,facies,'EdgeAlpha',.1);
caxis([.5 6.5]); colormap(gca,flipud(lines(6)));
mrstColorbar(facies,'west')
view(30,30); axis tight off; 
title('Facies','FontSize',10,'FontWeight','normal');

subplot(2,2,2,'FontSize',12) 
yyaxis left, bar(ncell), axis tight
yyaxis right, plot(cumsum(ncell),'LineWidth',1)
title('Cell count per layer','FontSize',10,'FontWeight','normal');

subplot(2,2,4,'FontSize',12) 
yyaxis left, bar(pvols), set(gca,'YScale','log','YTick',10.^(-4:1))
yyaxis right, plot(cumsum(pvols),'LineWidth',1), axis tight
title('Pore volume per layer','FontSize',10,'FontWeight','normal');

%% Set up and solve the fine-scale problem
% Pressure drop from left to right and no-flow conditions elsewhere,
% single-phase fluid model with unit density and viscosity
xf    = G.faces.centroids(:, 1);
left  = find(abs(xf-min(xf))< 1e-4);
right = find(abs(xf-max(xf))< 1e-4);
bc    = addBC([], left,  'pressure', 1);
bc    = addBC(bc, right, 'pressure', 0);

figure;
plotGrid(G, 'facecolor', 'none');
plotFaces(G, left, 'facecolor', 'r')
plotFaces(G, right, 'facecolor', 'b')
view(30, 30)

fluid = initSingleFluid('rho', 1, 'mu', 1);
hT    = computeTrans(G, rock);
state = initResSol(G, 0);
ref   = incompTPFA(state, G, hT, fluid, 'bc', bc);

% Extract fine-scale matrix to use later for computing basis functions
A  = getIncomp1PhMatrix(G, hT, state, fluid);

%% Compute multiscale solution(s)
% This model is challenging to partition because of the many eroded and
% semi-eroded layers, which we already have seen give large variations in
% cell sizes as well as in the number of cells and the rock volume in each
% layer. Many of the thin cells are strongly curved (and degenerate) so
% that the centroids fall outside the cell volume. Straightforward
% approaches, like a load-balanced partition in index space or using cell
% coordinates (centroids or vertices) to partition the cells into
% rectangular boxes in physical space can in some cases lead to grid blocks
% that give strong monotonicity violations. 
%
% We thereofre investigate an alternative approach in which we merge layers
% vertically so that each of the merged layers have a certain "thickness"
% measured in cell count or bulk/pore volume. The bar charts to the right
% in the first plot indicate that the layers consist of nine natural groups
% we can use as a starting point. Once we have obtained a satisfactory
% vertical partition, we can use a "cookie cutter" approach, in index or
% physical space, to partition in the horizontal direction.
n = 6;
[ms,part] = deal(cell(n,1));
for i=1:n
    switch i
        case 1    % 6x6x9 in physical space
            part{i} = sampleFromBox(G,reshape(1:6*6*9,[6,6,9]));
            part{i} = processPartition(G, compressPartition(part{i}));

        case 2    % 6x6x9 in index space
            part{i} = partitionUI(G,[6 6 9]);
            part{i} = processPartition(G, compressPartition(part{i}));

        case 3    % vertical by cell count + 6x6 horizontal
            cl = discretize(cumsum(ncell),linspace(0,G.cells.num+1,10));
            [~,edges] = rlencode(cl);
            part{i}   = partitionLayers(G, [6 6], [1; cumsum(edges)+1]);

        case 4    % vertical by pore volume + 6x6 horizontal
            cl = discretize(cumsum(pvols),linspace(0,sum(pvols)*1.01,10));
            [~,edges] = rlencode(cl);
            part{i}   = partitionLayers(G, [6 6], [1; cumsum(edges)+1]);
 
        case 5    % vertical by bulk volume + 6x6 horizontal
            cl = discretize(cumsum(vols),linspace(0,sum(vols)*1.01,10));
            [~,edges] = rlencode(cl);
            part{i}   = partitionLayers(G, [6 6], [1; cumsum(edges)+1]);

        case 6    % custom vertical + 6x6 horizontal
            L = [1 20 58 100 145 185 220 265 305 334]';
            part{i}   = partitionLayers(G, [6 6], L);
    end
    CG = generateCoarseGrid(G, part{i});
    CG = coarsenGeometry(CG);
    CG = storeInteractionRegion(CG,'edgeBoundaryCenters', true, 'adjustCenters', true);

    % Compute basis functions
    basis = getMultiscaleBasis(CG, A, 'type', 'msrsb', 'iterations', 150);

    % Compute multiscale solution
    ms{i} = incompMultiscale(state, CG, hT, fluid, basis, 'bc', bc);
end

%% Visualize two of the solutions
figure('Position',[300 200 1000 280])
titles = {'Partition in physical space', 'Partition in index space', ...
    'Partition by cell count', 'Partition by pore volume', ...
    'Partition by bulk volume', 'Manual partition'};
selection = [1 3];
for i=1:2
    subplot(1,3,i)
    plotCellData(G, ms{selection(i)}.pressure, 'EdgeAlpha',.1);
    view(30,30), axis tight off
    outlineCoarseGrid(G, part{selection(i)},'FaceColor','none');
    title(titles{selection(i)},'FontWeight','normal','FontSize',12)
    caxis([-1/32 33/32]);
end
colormap([1 1 1; parula(32); .7 .7 .7]);

% show the cell count partition
subplot(1,3,3) 
cl = discretize(cumsum(ncell),linspace(0,G.cells.num+1,10));
[~,edges] = rlencode(cl);
L = [1; cumsum(edges)+1];
area(cumsum(ncell));
hold on, 
for i=1:numel(L)
    plot(L([i i]),[1 G.cells.num],'r--','LineWidth',1);
end
axis tight, set(gca,'FontSize',12);

%% Compute discrepancies and monotonicity violations
pv      = poreVolume(G,rock);
errinf  = @(x) abs(x.pressure - ref.pressure)/(max(ref.pressure) - min(ref.pressure));
err2    = @(x) sum((x.pressure - ref.pressure).^2.*pv)./sum(pv);
labels = {'physical','index','cell count','pore vol','volume','manual'};
fprintf(1,'\n%-10s%-8s\t%-8s\t%-8s\t%-8s\n','method','L2','mean','max','std');
l2err   = zeros(n,1);
linf    = zeros(n,1);
violate = zeros(n,2);
range   = zeros(n,2);
for i=1:n
    l2err(i) = err2(ms{i});
    err      = errinf(ms{i});
    linf(i)  = max(err);
    violate(i,1) = sum(ms{i}.pressure<0);
    violate(i,2) = sum(ms{i}.pressure>1);
    range(i,1)   = min(ms{i}.pressure);
    range(i,2)   = max(ms{i}.pressure);
    fprintf(1,'%-10s%.2e\t%.2e\t%.2e\t%.2e\n', labels{i},...
        err2(ms{i}), mean(err), max(err), std(err));
end
figure('Position',[300 200 1000 340])
subplot(1,3,1)
semilogy(linf,'-o','LineWidth',1,'MarkerFaceColor',[.6 .6 .6]);
set(gca,'XTick',1:n,'XTickLabel',labels,'XTickLabelRotation',45, ...
    'FontSize',12,'Xlim',[.5 n+.5]);
title('L^\infty error', 'FontWeight','normal'), 

subplot(1,3,2)
semilogy(l2err,'-o','LineWidth',1,'MarkerFaceColor',[.6 .6 .6]);
set(gca,'XTick',1:n,'XTickLabel',labels,'XTickLabelRotation',45, ...
    'FontSize',12,'Xlim',[.5 n+.5]);
title('L^2 error', 'FontWeight','normal'), 

subplot(1,3,3)
barh(violate,'stacked');
set(gca,'YTick',1:n,'YTickLabel',labels, 'YDir', 'reverse',...
    'FontSize',12,'Ylim',[.25 n+.75]);
legend('<0','>1','location','southoutside','orientation','horizontal');
for i=1:n
    h=text(min(sum(violate(i,:))+10,220),i,...
        sprintf('[%.3f,%.3f]',range(i,1),range(i,2)));
end
title('# out-of-bound cells', 'FontSize', 12, 'FontWeight','normal');

%% Show the worst monotonicity violation
figure
sol   = ms{2}.pressure;
p     = part{2};
ind   = unique(p(sol>1.1));
neigh = getCellNeighbors(CG,ind(1));
plotGrid(CG,ind(2),'FaceAlpha',.2,'FaceColor',[.8 .3 .3]);
plotGrid(CG,ind(1),'FaceAlpha',.1,'FaceColor',[.3 .3 .8]);
plotCellData(G,sol,sol>1.1);
plotGrid(CG,neigh,'FaceAlpha',.1);
plotGrid(CG,'FaceAlpha',0,'EdgeAlpha',.3)
axis tight off;
axis([15 30 0 15 0 3]); view(-20,35)
title('Two blocks containing cells with highest discrepancy')

%% Copyright notice

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
