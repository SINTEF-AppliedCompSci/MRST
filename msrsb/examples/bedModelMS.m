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
% problems for the MsFV method. Here, we demonstrate that this is not the
% case for the MsRSB method.

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

subplot(2,2,3,'FontSize',12)
plotCellData(G,facies,'EdgeAlpha',.1);
caxis([.5 6.5]); colormap(gca,flipud(lines(6)));
mrstColorbar(facies,'west')
view(30,30); axis tight off;

subplot(2,2,2,'FontSize',12) 
yyaxis left, bar(ncell), axis tight
yyaxis right, plot(cumsum(ncell),'LineWidth',1)

subplot(2,2,4,'FontSize',12) 
yyaxis left, bar(pvols), set(gca,'YScale','log','YTick',10.^(-4:1))
yyaxis right, plot(cumsum(pvols),'LineWidth',1), axis tight

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
% layer. However, more important, many of the thin cells are strongly
% curved (and degenerate) so that the centroids fall outside the cell
% volume. Altogether, this means that straightforward approaches, like a
% load-balanced partition in index space or using cell coordinates
% (centroids or vertices) to partition the cells into rectangular boxes in
% physical space only work well for certain coarsening parameters. A better
% approach would be to start by merging layers vertically so that each of
% the merged layers have a certain "thickness" measured in cell count or
% bulk/pore volume. The bar charts to the right in the first plot indicate
% that the layers consist of nine natural groups we can use as a starting
% point. Once we have obtained a satisfactory vertical partition, we can
% use a "cookie cutter" approach, in index or physical space, to partition
% in the horizontal direction.
n = 6;
[ms,part] = deal(cell(n,1));
for i=1:n
    switch i
        case 1    % 6x6x5 in physical space
            part{i} = sampleFromBox(G,reshape(1:6*6*5,[6,6,5]));

        case 2    % 6x6x5 in index space
            part{i} = partitionUI(G,[6 6 5]);

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
sols = [{ref}; ms(:)];
part = [{ones(G.cells.num,1)}; part(:)];

%% Visualize the solutions
figure('Position',[300 200 1000 420])
pos = [1 2 3 5 6];
titles = {'fine scale', '6x6x5 physical space', '6x6x5 index space', ...
    'vertical by cell count', 'vertical by pore volume'};
for i=1:4
    subplot(2,3,pos(i),'FontSize',12)
    plotCellData(G, sols{i}.pressure, 'EdgeAlpha',.1);
    view(30,30), axis tight off
    if i>1, outlineCoarseGrid(G, part{i}); end
    title(titles{i},'FontWeight','normal')
end
subplot(2,3,4,'FontSize',12)
cl = discretize(cumsum(ncell),linspace(0,G.cells.num+1,10));
[~,edges] = rlencode(cl);
L = [1; cumsum(edges)+1];
bar(cumsum(ncell));
hold on, 
for i=1:numel(L)
    plot(L([i i]),[1 G.cells.num],'r--','LineWidth',1);
end
axis tight

%% Compute discrepancies
pv      = poreVolume(G,rock);
errinf  = @(x) abs(x.pressure - ref.pressure)/(max(ref.pressure) - min(ref.pressure));
err2    = @(x) sum((x.pressure - ref.pressure).^2.*pv)./sum(pv);
fprintf(1,'\n%-8s\t%-8s\t%-8s\t%-8s\n','L2','mean','max','std');
l2err = zeros(n,1);
for i=1:n
    l2err(i) = err2(ms{i});
    err      = errinf(ms{i});
    fprintf(1,'%.2e\t%.2e\t%.2e\t%.2e\n',...
        err2(ms{i}), mean(err), max(err), std(err));
end
subplot(2,3,6,'FontSize',12)
labels = {'physical','index','cells','pore vol','volume','manual'};
semilogy(l2err,'-o','LineWidth',1,'MarkerFaceColor',[.6 .6 .6]);
set(gca,'XTick',1:n,'XTickLabel',labels,'XTickLabelRotation',45);
title('L^2 error', 'FontWeight','normal')
%% Copyright notice

% <html>
% <p><font size="-1">
% Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.
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
