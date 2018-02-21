%% Generate Coarse Grids with Near-Well Refinement
% Pressure gradients and flow rates will typically be much larger near
% wells than inside the reservoir. The accuracy with which we represent the
% flow in and out of wells will to a large extent determine the accuracy of
% an overall simulation and as a result one therefore often desires to have
% higher grid resolution in the near-well zone than inside the reservoir.
%
% The example considers the smallest pillar grid from ModelB4, which can
% either be partitioned rectangularly in index space, or using METIS with
% transmissibilities as edge weights. The latter approach gives blocks with
% boundaries aligning to sharp media contrasts. On top of this, we split
% partitions across faults and use the function 'refineNearWell' to impose
% radial refinement in near-well regions.

useMetis = false;

%% Load grid model and create petrophysical data
file = fullfile(getDatasetPath('CaseB4'),'stairstep_36x48.grdecl');
G = processGRDECL(readGRDECL(file));
G = computeGeometry(G);
layers = [1 2 7 12 13];
rock.perm = logNormLayers(G.cartDims, [50 300 100 20]*milli*darcy, ...
    'sz', [51 3 3], 'indices', layers);

W = verticalWell([], G, rock, G.cartDims(1), 1, [], 'InnerProduct', 'ip_tpf');
W = verticalWell(W, G, rock, 14, 36, [],'InnerProduct', 'ip_tpf');

clf
plotCellData(G, log10(rock.perm),'EdgeAlpha',.1); view(3); 
plotWell(G, W, 'height', 30,'FontSize',12); 
colormap(parula), axis tight off, view(250, 50)

%% Uniform/METIS partition + fault faces
% We consider two possible partitions. The first is a standard
% load-balanced partition in index space. The second is an unstructured
% partitioning by the METIS graph library that adapts to the underlying
% geology. To this end, we transmissibilities as edge-weights in the
% graph-partitioning algorithm of METIS so that it tries to make grid
% blocks having as homogeneous permeability as possible. For both
% partitions, we also split blocks across faults.
if useMetis
    T = computeTrans(G, rock);
    p0 = partitionMETIS(G, T, 7*9*4, 'useLog', true);
else
    p0 = partitionUI(G, [7 9 1]); max(p0)
end
p0 = processPartition(G, p0, find(G.faces.tag==1)); max(p0)
plotPartition = @(p) plotCellData(G, p, 'EdgeColor','k','EdgeAlpha',.05);
clf
plotPartition(p0); outlineCoarseGrid(G, p0)
plotWell(G, W, 'height', 30,'FontSize',12)
axis tight off, view(250, 50),
cmap = tatarizeMap(max(p0)); colormap(.5*cmap+.5*ones(size(cmap)));

%% Radial refinement
% The function 'refineNearWell' takes a set of points and partitions these
% according to the distance in the xy-plane from a single well point. Here,
% we will this function to refine the coarse blocks that contains wells.
% The first well is placed in the corner and we partition the perforated
% well block into five radial sections. The width of the radial sections is
% set to decay as log(r). For the second well, we refine all the
% neighboring blocks surrounding the well block using a number of angular
% sectors that increases as we move radially out from the well point.
CG  = generateCoarseGrid(G, p0);
clf, plotGrid(CG); view(3);
angSectors = {1, [1,4,6,9]};
radSectors = [4 4];
p = p0;
for i=1:numel(W)
    wc    = W(i).cells(1);
    wpt   = G.cells.centroids(wc,:);
    pw    = p(wc);
    if i==1
        % Pick all cells inside the well block
        cells = (p==pw);
    else
        % Use adjacency matrix to compute nearest neighbors. We start with
        % a vector e which is equal one in the well block and zero
        % elsewhere. Multiplying by the adjacency matrix A will set the
        % value in each block equal the sum over the block and its
        % face-neighbors. After two multiplications, all blocks surrounding
        % the initial block should have value larger than unity.
        A     = getConnectivityMatrix(getNeighbourship(CG),true,CG.cells.num);
        rblk  = zeros(CG.cells.num,1); rblk(pw)=1; rblk((A*A*rblk)>1) = 1;
        cells = rblk(p)>0;
    end
    pts   = G.cells.centroids(cells,:);
    out   = refineNearWell(pts, wpt, 'angleBins', angSectors{i}, ...
        'radiusBins', radSectors(i), 'logbins', true, 'maxRadius', inf);
    p(cells) = max(p) + out;
end
p = compressPartition(p);
clf
plotPartition(p); outlineCoarseGrid(G, p)
plotWell(G, W, 'height', 30,'FontSize',12)
axis tight off, view(250, 50), 
cmap = tatarizeMap(max(p)); colormap(.5*cmap+.5*ones(size(cmap)));

%% Coarsen in the vertical direction
% For the Cartesian partition, we also coarsen in the vertical direction
% following the layering in the petrophysical data, which we assume is
% known.
if ~useMetis
    pK = rldecode((1:numel(layers)-1)',diff(layers));
    [~,~,k]=gridLogicalIndices(G);
    p = compressPartition((pK(k)-1)*max(p)+p);
    p = processPartition(G, p);
    rp = randperm(max(p))';
    clf,
    plotPartition(rp(p)); outlineCoarseGrid(G, p)
    plotWell(G, W, 'height', 30,'FontSize',12)
    axis tight off, view(250, 50),
    cmap = tatarizeMap(max(p)); colormap(.5*cmap+.5*ones(size(cmap)));
end