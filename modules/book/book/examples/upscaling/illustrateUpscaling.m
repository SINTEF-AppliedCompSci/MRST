%% Figures for conceptual illustration of upscaling
% We generate a grid that consists of two facies and then perform an
% upscaling so that 5x5x5 fine cells are replaced by a coarse block. We
% show the following four plots: fine grid with one coarse cell,
% porosity inside this coarse cell, the porosity inside the coarse block,
% and the coarse model

mrstModule add coarsegrid upscaling

%% Fine model
G  = computeGeometry(cartGrid([40 20 15]));
K1 = gaussianField(G.cartDims, [300 2000]); 
p1 = K1(:)*1e-4 + .2;
K2 = gaussianField(G.cartDims, [10 900]);
p2 = K2(:)*1e-4 + .2;

rad1 = G.cells.centroids(:,1).^2 + .5*G.cells.centroids(:,2).^2 ...
   + (G.cells.centroids(:,3)-2).^2;
rad2 = .5*(G.cells.centroids(:,1)-40).^2 + 2*G.cells.centroids(:,2).^2 ...
   + 2*(G.cells.centroids(:,3)-2).^2;

ind = ((rad1>600) & (rad1<1500)) | ((rad2>700) & (rad2<1400));
rock.poro = p2(:);
rock.poro(ind) = p1(ind);

% Plot fine model
figure('Position', [680 280 560 540]);
plotCellData(G,rock.poro,'EdgeColor','k','edgealpha',.1); 
view(3); axis tight off
[hc,hh]=colorbarHist(rock.poro,[.2 .4],'South', 80);
pos=get(hc,'Position'); set(hc,'Position',pos - [.02 0 .2 -.01],'FontSize',12);
pos=get(hh,'Position'); set(hh,'Position',pos - [.02 0.02 .2 -.01]);
set(gca,'Position',[.13 .2 .775 .775])

% Generate and plot coarse grid block
p = partitionCartGrid(G.cartDims, G.cartDims/5);
CG = generateCoarseGrid(G, p);
crock.poro = accumarray(p, rock.poro)./accumarray(p,1);
plotGrid(CG,1,'FaceColor','none','EdgeColor','k','LineWidth',1.5);
set(gcf,'PaperPositionMode','auto');
% print -dpng illUpscaling-fine.png;

%% Plot cells inside a single block
figure;
plotCellData(G,rock.poro, p==1, 'EdgeColor','k'); 
view(3); axis equal off; caxis([0.2 0.4]);
% print -dpng illUpscaling-cells.png;

%% Plot the corresponding coarse block
figure;
cG = computeGeometry(cartGrid(G.cartDims/5,G.cartDims));
plotCellData(cG,crock.poro,1,'EdgeColor','k','LineWidth',1.5);
view(3); axis equal off; caxis([0.2 0.4]);
% print -dpng illUpscaling-block.png;

%% Coarse model
figure('Position', [680 280 560 540]);
plotCellData(cG,crock.poro,'EdgeColor','k','edgealpha',.1); 
view(3); axis tight off

% Generate and plot coarse grid block
plotGrid(CG,1,'FaceColor','none','EdgeColor','k','LineWidth',1.5);

[hc,hh]=colorbarHist(crock.poro,[.2 .4],'South', 40);
pos=get(hc,'Position'); set(hc,'Position',pos + [.2    0 -.2 .01],'FontSize',12);
pos=get(hh,'Position'); set(hh,'Position',pos + [.2 -.025 -.2 .01]);
set(gca,'Position',[.13 .2 .775 .775])
set(gcf,'PaperPositionMode','auto');
% print -dpng illUpscaling-coarse.png;