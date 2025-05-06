%% Exercise 3.5.1
% Populate grids from Figure 3.31 with petrophysical data that preserves
% the layering. These are corner-point grids and are simple to populate
grdecl = simpleGrdecl([20, 20, 5],@(x) .25*(x-.5), 'flat' ,true);
G = computeGeometry(processGRDECL(grdecl));
K = logNormLayers(G.cartDims,100*randperm(5)*milli*darcy,'sigma',1.5);
clf, plotCellData(G,log10(K)); view(3);

%%
% Or using a subsample from the SPE10 data set
mrstModule add spe10
rock = getSPE10rock(1:20,1:20,33:37);
clf, plotCellData(G,rock.poro); view(3);

%%
% Likewise for the grid shown to the left in Figure 3.32
grdecl=makeModel3([30,20,5]);
G = computeGeometry(processGRDECL(grdecl));
K = logNormLayers(G.cartDims,100*randperm(5)*milli*darcy,'sigma',1.5);
clf, plotCellData(G,log10(K(G.cells.indexMap))); view(50,30);

%%
% The grids shown in the middle and right plots of Figure 3.32 are not
% necessarily easy to populate because they do not contain explicit layer
% information.
bdim = [25 25 5];
G = computeGeometry(extrudedTriangleGrid(50,true));
Kx = reshape(...
    logNormLayers(bdim,50*(1:bdim(3))*milli*darcy,'sigma',0.5), bdim);
p = sampleFromBox(G,Kx);
clf, plotCellData(G, p), view(-160,20);

%%
% We can make an attempt to reconstruct the layering from the information
% that is available in the grid structure, i.e., using G.numLayers and
% G.layerSize
layerID = repmat((1:G.numLayers)',1,G.layerSize)';
clf, plotCellData(G, layerID(G.cells.indexMap)), view(-160,20);

%%
% If we now temporarily substitute layerID as the z-component of the cell
% centroids, we can once again try to sample from the bounding box
pts = G.cells.centroids;
pts(:,3) = layerID(G.cells.indexMap);
[~,i]  = min(pts,[],1);
[~,j]  = max(pts,[],1);
n      = unique(gridCellNodes(G,[i j]));
Y      = G.nodes.coords(n,:);
mY     = min(Y); mY(3) = .5;
MY     = max(Y); MY(3) = G.numLayers+.5;
coords = bsxfun(@minus, pts, mY);
coords = bsxfun(@rdivide, coords, MY - mY);
ijk    = num2cell(ceil(bsxfun(@times, coords, bdim)),1);
p      = Kx(sub2ind(bdim,ijk{:}));
clf, plotCellData(G, p), view(-160,20); axis tight

%%
% We can also try to reconstruct the correct layerID's in the case that
% there has been displacement along the fault so that the top cell layers
% on both sides represent the same geological layer
arealID = repmat((1:G.layerSize)',G.numLayers,1);
ind = zeros(G.numLayers,G.layerSize);
ind(sub2ind(size(ind),layerID(G.cells.indexMap), arealID(G.cells.indexMap))) = 1;
ind = cumsum(ind).';
clf, plotCellData(G, ind(G.cells.indexMap)); view(-160,20); axis tight


%% Exercise 3.5.2
pth = getDatasetPath('BedModel2');
grdecl = readGRDECL(fullfile(pth,'BedModel2.grdecl'));
grdecl = convertInputUnits(grdecl, getUnitSystem('METRIC'));

% First, show active and inactive cells
Gc = cartGrid(grdecl.cartDims);
clf
plotCellData(Gc,grdecl.ACTNUM,'EdgeAlpha',.05); view(3); axis tight
hc=colorbar('horiz'); colormap(winter(2)); caxis([-.5 1.5]);
set(hc,'XTick',0:1,'XTickLabel',{'Inactive','Active'});

%%
G      = computeGeometry(processGRDECL(grdecl));
rock   = grdecl2Rock(grdecl, G);
nfaces = diff(G.cells.facePos);
clf
plotGrid(cartGrid(ones(1,3), grdecl.cartDims),'FaceColor','none');
plotCellData(Gc, nfaces, G.cells.indexMap, 'EdgeAlpha',.05); 
view(3), axis tight,
mC = min(nfaces); MC = max(nfaces); colormap(jet(MC-mC+1));
caxis([mC-.5 MC+.5]); hc=colorbar; set(hc,'YTick',mC:MC)

%%
clf
plotCellData(G, nfaces,'Edgealpha',.05); view(3), axis tight; colorbar
caxis([mC-.5 MC+.5]); hc=colorbar; set(hc,'YTick',mC:MC)


%% Exercise 3.5.3
% Not available here


%% Exercise 3.5.4
% Not available here


%% Exercise 3.5.5
% Not available here