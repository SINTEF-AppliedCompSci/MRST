%%
mrstModule add agglom coarsegrid diagnostics incomp spe10

%% CaseB4: merge small cells
% Load and set up model
file = fullfile(getDatasetPath('CaseB4'),'pillar_36x48.grdecl');
G = processGRDECL(readGRDECL(file));
G = computeGeometry(G);
layers = [1 2 8 12 13];
K = logNormLayers(G.cartDims, [50 300 100 20]*milli*darcy, ...
                  'sz', [51 3 3], 'indices', layers);
rock = makeRock(G, bsxfun(@times, K*milli*darcy,[1 1 .1]), 0.2);

W = verticalWell([], G, rock, G.cartDims(1), 1, [], 'InnerProduct', 'ip_tpf');
W = verticalWell(W, G, rock, 14, 36, [],'InnerProduct', 'ip_tpf');

% Partition topology and split across faults 
p0 = partitionUI(G, [7 9 1]);
pf = processPartition(G, p0, find(G.faces.tag==1));
figure, myPlotPartition(G, W, pf, 400);

%%
% Merge small blocks
% To avoid merging across faults, we generate a static partition that
% splits the model into three different fault blocks
I  = poreVolume(G, rock)./G.cells.volumes;
pm = mergeBlocks(pf, G, I, I, 50, 'static_partition', ...
    processPartition(G,ones(G.cells.num,1),find(G.faces.tag==1)));
figure, myPlotPartition(G, W, pm, 400);

%% Adaptive Cartesian coarsening for layer of SPE10
% We consider a subset of layer 25 from SPE 10 and study how recursive
% Cartesian refinement works for three different indicator functions:
% permeability, velocity, and time of flight.

% Setup model
[G, W, rock] = getSPE10setup(25);
for i = 1:numel(W), W(i).compi = 1; end
rock.poro = max(rock.poro, 1e-4);
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);
hT = computeTrans(G, rock);
rS = incompTPFA(initState(G,W,0), G, hT, fluid, 'wells', W);
Tf = computeTimeOfFlight(rS, G, rock, 'wells', W);
Tb = computeTimeOfFlight(rS, G, rock, 'wells', W, 'reverse', true);
pv = poreVolume(G, rock);

% Construct indicators
iK = log10(rock.perm(:,1)); iK = iK - min(iK) + 1;
v  = sqrt(sum(faceFlux2cellVelocity(G, rS.flux).^2, 2));
iV = log10(v); iV = iV - min(iV) + 1;
iT = -log10(Tf+Tb); iT = iT - min(iT) + 1;

% Extract subset of model
p1  = partitionUI(G, [3, 11, 1]);
ind = find(p1<19);
G   = extractSubgrid(G, ind);
[p1,iK,iV,iT,pv] = deal(p1(ind),iK(ind),iV(ind),iT(ind),pv(ind));

% Recursively subdivide the blocks by a factor 2x2
pK = refineUniformShape(p1, G, iK, 25, 'CartDims', [2,2,1]);
pV = refineUniformShape(p1, G, iV, 25, 'CartDims', [2,2,1]);
pT = refineUniformShape(p1, G, iT, 25, 'CartDims', [2,2,1]);

% Plot the results
I = {iK, iV, iT}; p = {pK, pV, pT}; 
head = {'Permeability', 'Velocity', 'Time-of-flight'};
wc = vertcat(W.cells);
x = G.cells.centroids(wc(wc<G.cells.num),:);
clf; set(gcf,'Position',[510 420 850 400]);
for i=1:3
    subplot(1,3,i)
    plotCellData(G,I{i},'EdgeColor','none');
    hold on
    plot(x(:,1),x(:,2),'ok','MarkerSize',7,'MarkerFaceColor','w'); hold off
    outlineCoarseGrid(G, p{i}, 'EdgeColor','k'); axis tight off
    title([head{i} ' indicator']);
end
colormap(jet(128)*.7+.3*ones(128,3));


%% Non-uniform coarsening
% We outline the main steps of the original NUC algorithm proposed by
% Aarnes et al. To this end, we extract a further subdivision of the model
ind = find(p1>9);
G   = extractSubgrid(G, ind);
[p1,iK,iV,iT,pv] = deal(p1(ind),iK(ind),iV(ind),iT(ind),pv(ind));
volI = pv./G.cells.volumes;
flwI = iV;

%%
clf,set(gcf,'Position',[50 560 1400 240]);
subplot(1,6,1)
plotCellData(G, flwI, 'EdgeColor','none'); axis tight
set(gca,'XTick',[],'YTick',[]); title('Indicator');

subplot(1,6,2)
ps = segmentIndicator(G, flwI, 5);
plotCellData(G, ps, 'EdgeColor','none'); axis tight
outlineCoarseGrid(G, ps, 'EdgeColor','k','FaceColor','none');
set(gca,'XTick',[],'YTick',[]); title('segmentIndicator');
xlabel(['#blocks: ' num2str(max(ps))]);

subplot(1,6,3)
pm1 = mergeBlocks(ps, G, volI, flwI, 20);
plotCellData(G, pm1, 'EdgeColor','none'); axis tight
outlineCoarseGrid(G, pm1, 'EdgeColor','k','FaceColor','none');
set(gca,'XTick',[],'YTick',[]); title('mergeBlocks');
xlabel(['#blocks: ' num2str(max(pm1))]);

subplot(1,6,4)
pr = refineGreedy(pm1, G, flwI, 30);
plotCellData(G, pr, 'EdgeColor','none'); axis tight
outlineCoarseGrid(G, pr, 'EdgeColor','k','FaceColor','none');
set(gca,'XTick',[],'YTick',[]); title('refineGreedy');
xlabel(['#blocks: ' num2str(max(pr))]);

subplot(1,6,5)
pm2 = mergeBlocks(pr, G, volI, flwI, 20);
plotCellData(G, pm2, 'EdgeColor','none'); axis tight
outlineCoarseGrid(G, pm2, 'EdgeColor','k','FaceColor','none');
set(gca,'XTick',[],'YTick',[]); title('mergeBlocks');
xlabel(['#blocks: ' num2str(max(pm2))]);

subplot(1,6,6);
p = mergeBlocks2(ps, G, volI, flwI, 20, 30);
p = refineGreedy2(p, G, flwI, 30, 'nlevel',1);
p = mergeBlocks2(p, G, volI, flwI, 20, 30);
plotCellData(G, p, 'EdgeColor','none'); axis tight
outlineCoarseGrid(G, p, 'EdgeColor','k','FaceColor','none');
set(gca,'XTick',[],'YTick',[]); title('Version 2');
xlabel(['#blocks: ' num2str(max(p))]);

colormap(jet(128)*.6+.4*ones(128,3));
