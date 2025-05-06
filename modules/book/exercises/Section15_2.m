%% Exercise 15.2.1
% Upscale additive properties from SAIGUP
mrstModule add libgeometry deckformat coarsegrid
grdecl=readGRDECL(fullfile(getDatasetPath('SAIGUP'),'SAIGUP.GRDECL'));
G = mcomputeGeometry(processgrid(grdecl));
rock = grdecl2Rock(grdecl,G.cells.indexMap);
p = partitionUI(G, G.cartDims./[5 5 4]);
p = compressPartition(processPartition(G,p));
crock1.poro = accumarray(p,rock.poro.*G.cells.volumes)./...
   max(accumarray(p,G.cells.volumes),eps);
crock1.ntg = accumarray(p,rock.poro.*G.cells.volumes.*rock.ntg)./...
   accumarray(p,G.cells.volumes.*rock.poro);

%%
% check that we have upscaled correctly
pv = poreVolume(G,rock);
CG1 = coarsenGeometry(generateCoarseGrid(G,p));
max(abs(accumarray(CG1.partition,poreVolume(G,rock))-poreVolume(CG1,crock1)))

%%
% Visualize the fine and coarse model
figure
subplot(2,2,1);
plotCellData(G,rock.poro,'EdgeColor','k','EdgeAlpha',.1);
view(-98,50), axis tight off; cx = caxis();
subplot(2,2,2);
plotCellData(CG1,crock1.poro,'EdgeColor','none');
plotGrid(CG1,'EdgeColor','k','FaceColor','none','EdgeAlpha',.1);
view(-98,50), axis tight off; caxis(cx);

%%
% Make grid that adapts to SATNUM
sn = grdecl.SATNUM(G.cells.indexMap);
p = partitionUI(G,[G.cartDims(1:2)/5 1]);
p = compressPartition(processPartition(G, p + max(p)*sn));
CG2 = coarsenGeometry(generateCoarseGrid(G,p));
crock2.poro = accumarray(p,rock.poro.*G.cells.volumes)./...
   max(accumarray(p,G.cells.volumes),eps);

subplot(2,2,3); cla
plotCellData(CG2,crock2.poro,'EdgeColor','none');
plotGrid(CG2,'EdgeColor','k','FaceColor','none','EdgeAlpha',.1);
view(-98,50), axis tight off; caxis(cx);

%%
% Make grid with original vertical resolution
p = partitionUI(G, G.cartDims./[5 5 1]);
p = compressPartition(processPartition(G,p));
crock3.poro = accumarray(p,rock.poro.*G.cells.volumes)./...
   max(accumarray(p,G.cells.volumes),eps);
crock3.ntg = accumarray(p,rock.poro.*G.cells.volumes.*rock.ntg)./...
   accumarray(p,G.cells.volumes.*rock.poro);
CG3 = coarsenGeometry(generateCoarseGrid(G,p));

subplot(2,2,4);
plotCellData(CG3,crock3.poro,'EdgeColor','none');
plotGrid(CG3,'EdgeColor','k','FaceColor','none','EdgeAlpha',.1);
view(-98,50), axis tight off; caxis(cx);

%%
% show histograms
px = [min(rock.poro),max(rock.poro)];
bins = linspace(px(1),px(2),51);
pargs = {'EdgeColor',[0.1 0.1 0.1],'FaceColor',[.6 .6 .6]};
figure
subplot(2,2,1); hist(rock.poro,bins);
set(get(gca,'Children'), pargs{:}); axis tight, set(gca,'XLim',px);
title('Fine scale');

subplot(2,2,2); hist(crock1.poro,bins);
set(get(gca,'Children'), pargs{:}); axis tight, set(gca,'XLim',px);
err1 = sum(abs(rock.poro - crock1.poro(CG1.partition)))/sum(rock.poro)*100;
title(['8x24x5: ' num2str(err1) '%'])

subplot(2,2,3); hist(crock2.poro,bins);
set(get(gca,'Children'), pargs{:}); axis tight, set(gca,'XLim',px);
err2 = sum(abs(rock.poro - crock2.poro(CG2.partition)))/sum(rock.poro)*100;
title(['8x24 + satnum: ' num2str(err2) '%'])

subplot(2,2,4); hist(crock3.poro,bins);
set(get(gca,'Children'), pargs{:}); axis tight, set(gca,'XLim',px);
err3 = sum(abs(rock.poro - crock3.poro(CG3.partition)))/sum(rock.poro)*100;
title(['8x24x20: ' num2str(err3) '%']);

%% Exercise 15.2.2
% Upscale model from BedModels1 dataset
grdecl=readGRDECL(fullfile(getDatasetPath('BedModels1'),'mortarTestModel.grdecl'));
G = computeGeometry(processGRDECL(grdecl));
rock = grdecl2Rock(grdecl,G.cells.indexMap);

%%
% Visualize porosity
figure('Position', [680 280 560 540]);
plotCellData(G,rock.poro,'EdgeColor','k','edgealpha',.1); 
view(3); axis tight off
set(gca,'Position',[.13 .2 .775 .775])
pax = [min(rock.poro) max(rock.poro)];
colorbar('SouthOutside','Position',[.1 .05 .6 .03],'FontSize',12);
caxis(pax);
axes('Position',[.1 .09 .6 .1]);
hist(rock.poro,linspace(pax(1),pax(2),101)); axis tight off
set(gca,'XLim',pax)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','none','EdgeColor',[.1 .1 .1]); 

%% 
% Make coarse grid
p = partitionUI(G,[5 5 1]);
q = partitionUI(G,[1 1 23]);
lavg = accumarray(q,rock.poro)./accumarray(q,1);
lfac = zeros(size(lavg));
lfac(lavg>.05) = 1;
lfac(lavg>.25) = 2;
lfac = compressPartition([0; cumsum(abs(diff(lfac)))]);
p = processPartition(G, p + max(p)*lfac(q));

CG = coarsenGeometry(generateCoarseGrid(G, p));
figure; 
plotCellData(G,p,'EdgeColor','none');
plotGrid(CG,'FaceColor','none','EdgeColor','k');
view(3);

%%
% average and visualize result
crock.poro = accumarray(p,rock.poro.*G.cells.volumes)./ ...
   accumarray(p,G.cells.volumes);
figure,
plotCellData(CG, crock.poro,'EdgeColor','none');
plotGrid(CG,'FaceColor','none','EdgeColor','k');
view(3);

%%
% Illustrate problems with flow simulation
figure,
I = 451:CG.cells.num;
I = I(CG.cells.centroids(I,3)>-.01);
plotCellData(CG,crock.poro(I),I);
plotGrid(CG, I,'EdgeColor','k','FaceColor','none');
view(3);

