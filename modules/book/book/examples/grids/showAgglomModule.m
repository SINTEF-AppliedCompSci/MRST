%%
mrstModule add agglom coarsegrid

%% ------ S E G M E N T A T I O N ------ 
% Demonstrate segmentation: apply to CaseB4 model
rng(1,'twister');    % reset random generator to reproduce the same results
file = fullfile(getDatasetPath('CaseB4'),'pillar_36x48.grdecl');
G = processGRDECL(readGRDECL(file));
G = computeGeometry(G);
layers = [1 2 8 12 13];
K = logNormLayers(G.cartDims, [50 300 150 20],'sz', [51 3 3], 'indices', layers);

%% 
% Plot permeability and outline initial bins
% We segment the permeability, but only out put the initial bins without
% splitting disconnected components
p = segmentIndicator(G, K, [0 30 80 205 inf],'split',false);

figure(1); clf
plotCellData(G, K,'EdgeAlpha',.1);
view(-120,15), axis tight, set(gca,'Projection','Perspective');
colormap(.8*jet+0.2*ones(size(jet)));
[~, hh]=colorbarHist(K, [0 450], 'South', 100);
h = outlineCoarseGrid(G,p,'EdgeColor','k','LineWidth',2,'FaceColor','none');
hold(hh,'on');
ym = .8*max(reshape(get(get(hh,'Children'),'YData'),[],1));
plot(hh,[30 80 205; 30 80 205], [0 0 0; ym ym ym],'Color',[.3 .3 .3],'LineWidth',2);
hold(hh,'off'); axis off
%set(hc,'FontSize',16);

%%
% Plot bins resulting from the segmentation
figure(2); clf
plotCellData(G, p,'EdgeAlpha',.1);
mp = max(p);
caxis([.5 mp+.5]); colormap(.8*tatarizeMap(mp)+.2*ones(mp,3));
view(-120,15), axis tight off, set(gca,'Projection','Perspective');
[hc,hh] = colorbarHist(p, [.5 mp+.5], 'South', mp);
bar(hh,accumarray(p,1),'FaceColor','none');
set(hh,'XLim',[.5 mp+.5]);
axis(hh,'off');
set(hc,'XTick',1:4); % set(hc,'FontSize',16);

%% ----- S A N I T Y   C H E C K S ------
% Split connected components (done by initial segmentation)
p = segmentIndicator(G, K, [0 30 80 205 inf]);
figure(1); clf
plotCellData(G, p,'EdgeAlpha',.1);
mp = max(p);
caxis([.5 mp+.5]); colormap(.8*tatarizeMap(mp)+.2*ones(mp,3));
view(-120,15), axis tight off, set(gca,'Projection','Perspective');
[hc,hh] = colorbarHist(p, [.5 mp+.5], 'South', mp);
bar(hh,log10(accumarray(p,1)),'FaceColor','none');
set(hh,'XLim',[.5 mp+.5]);
axis(hh,'off');
set(hc,'FontSize',16);

%%
% Detect confined blocks
cb = findConfinedBlocks(G, p);
flag=false(mp,1); flag(cb)=true;
plotGrid(G,flag(p),'FaceColor','k','EdgeColor','k');
delete([hc,hh]);

%% 
% Remove confined blocks
%[cb,p] = findConfinedBlocks(G,p);
figure(1); clf,
hp = plotCellData(G, p,'EdgeAlpha',.1);
mp = max(p);
caxis([.5 mp+.5]); colormap(.8*tatarizeMap(mp)+.2*ones(mp,3));
view(-120,15), axis tight off, set(gca,'Projection','Perspective');
[hc,hh] = colorbarHist(p, [.5 mp+.5], 'South', mp);
bar(hh,log10(accumarray(p,1)),'FaceColor','none');
set(hh,'XLim',[.5 mp+.5]);
axis(hh,'off');
set(hc,'FontSize',14);

%% Remove recursively confined blocks
% The implementation in findConfinedBlocks is relatively simple and does
% not work properly for recursivley confined blocks. For this, we can
% instead use graph algorithms from the MATLAB Boost Graph Library.
mrstModule add matlab_bgl
figure('Position',[450 540 940 260]);
G = computeGeometry(cartGrid([7 7]));
% p = ones(7,7); p(2:6,2:6)=2; p(3:5,3:5)=3;p(4,4)=4;
p = ones(7,7);p(:,4:7)=2;p(2:6,2:6)=3;p(3:5,3)=4;p(3:5,4)=5;p(3:5,5)=6;
pn = removeConfinedBlocks(G,p(:));
[~,pm] = findConfinedBlocks(G, p(:));
subplot(1,3,1), plotCellData(G,p(:)); axis equal off, caxis([1 7]), title('initial');
subplot(1,3,2), plotCellData(G,pm); axis equal off, caxis([1 7]), title('findConfinedBlocks');
subplot(1,3,3), plotCellData(G,pn);  axis equal off, caxis([1 7]), title('removeConfinedBlocks');

%% ------ M E R G I N G   F U N C T I O N S ------
G = computeGeometry(cartGrid([5 4]));
I = [1 1 5 2 2; 1 5 1 1 1; 4 5 1 2 2; 6 1 3 3 3]'; I=I(:);
p = segmentIndicator(G, I, 6);
n = G.cells.num;

set(figure(1),'Position',[450 540 940 260]);

subplot(1,3,1), cla
plotCellData(G,I,'edgecolor','none');
outlineCoarseGrid(G, p, 'k', 'LineWidth',2); axis off
colormap(.6*tatarizeMap(6)+.4*ones(6,3));
hold on
for i=1:n,
    text(G.cells.centroids(i,1),G.cells.centroids(i,2), num2str(I(i)));
end
hold off

% mergeBlocks
p1 = mergeBlocks(p, G, ones(n,1), I, 2);
subplot(1,3,2), cla
plotCellData(G,I,'EdgeColor','none');
outlineCoarseGrid(G, p1, 'k', 'LineWidth', 2); axis off
xm = sparse(p1,1:n,1)*[G.cells.centroids, ones(n,1), I];
Ib = xm(:,end);
xm = bsxfun(@rdivide, xm(:,1:2), xm(:,3));
for i=1:size(xm,1)
    text(xm(i,1),xm(i,2),num2str(Ib(i,end)));
end

% mergeBlocks2
p2 = mergeBlocks2(p, G, ones(G.cells.num,1), I, 2, 3);
subplot(1,3,3), cla
plotCellData(G,I,'EdgeColor','none');
outlineCoarseGrid(G, p2, 'k', 'LineWidth', 2); axis off
xm = sparse(p2,1:n,1)*[G.cells.centroids, ones(n,1), I];
Ib = xm(:,end);
xm = bsxfun(@rdivide, xm(:,1:2), xm(:,3));
for i=1:size(xm,1)
    text(xm(i,1),xm(i,2),num2str(Ib(i,end)));
end

%% ------ R E F I N E M E N T   M E T H O D S -----
pargs={'EdgeColor','none'};
blk = [5 13];
dim = {[4,4],[5,5]};
for n=1:2
    G = computeGeometry(cartGrid(dim{n}));
    I = ones(G.cells.num,1);

    figure('Position', [440 560 930 240]);

    subplot(1,6,1);
    plotGrid(G,'FaceColor','none'); axis equal off
    for i=1:G.cells.num,
        text(G.cells.centroids(i,1),G.cells.centroids(i,2),num2str(i),...
            'HorizontalAlignment','center','FontSize',8);
    end

    p = refineUniform(I,G,I,4,'CartDims',[2,2]);
    subplot(1,6,2); plotCellData(G,p,pargs{:}); axis equal off
    outlineCoarseGrid(G,p,'k'); caxis([.5 blk(n)+.5]);
    title('refineUniform','Position',[dim{n} 1].*[.5 1.02 1.02]);
        
    p = refineGreedy(I,G,I,4,'nlevel',n);
    subplot(1,6,3), plotCellData(G,p,pargs{:}); axis equal off,
    outlineCoarseGrid(G,p,'k');caxis([.5 blk(n)+.5]);
    title('refineGreedy','Position',[dim{n} 1].*[.5 1.02 1.02]);

    p = refineGreedy2(I,G,I,4,'nlevel',n);
    subplot(1,6,4), plotCellData(G,p,pargs{:}); axis equal off,
    outlineCoarseGrid(G,p,'k'); caxis([.5 blk(n)+.5]);
    title('refineGreedy2','Position',[dim{n} 1].*[.5 1.02 1.02]);

    p = refineGreedy3(I,G,I,4,'nlevel',n);
    subplot(1,6,5), plotCellData(G,p,pargs{:}); axis equal off,
    outlineCoarseGrid(G,p,'k');caxis([.5 blk(n)+.5]);
    title('refineGreedy3','Position',[dim{n} 1].*[.5 1.02 1.02]);

    p = refineGreedy4(I,G,I,4,'nlevel',n);
    subplot(1,6,6), plotCellData(G,p,pargs{:}); title('refineGreedy4');axis equal off,
    outlineCoarseGrid(G,p,'k');caxis([.5 blk(n)+.5]);
    title('refineGreedy4','Position',[dim{n} 1].*[.5 1.02 1.02]);
    
    colormap(.8*tatarizeMap(blk(n))+.2*ones(blk(n),3));
    h = colorbar('horiz'); 
    set(h,'Position',[0.29 0.13 0.455 0.1],'XTick',1:blk(n));
end
