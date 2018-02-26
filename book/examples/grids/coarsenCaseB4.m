%% Generate Coarse Grids with Near-Well Refinement
% Pressure gradients and flow rates will typically be much larger near
% wells than inside the reservoir. The accuracy with which we represent the
% flow in and out of wells will to a large extent determine the accuracy of
% an overall simulation and as a result one therefore often desires to have
% higher grid resolution in the near-well zone than inside the reservoir.
%
% The example considers the smallest pillar grid from CaseB4, which can
% either be partitioned rectangularly in index space, or using METIS with
% transmissibilities as edge weights. The latter approach gives blocks with
% boundaries aligning to sharp media contrasts. On top of this, we split
% partitions across faults and use the function 'refineNearWell' to impose
% radial refinement in near-well regions.
mrstModule add coarsegrid incomp
useMetis = false;

%% Load grid model and create petrophysical data
file = fullfile(getDatasetPath('CaseB4'),'pillar_36x48.grdecl');
G = processGRDECL(readGRDECL(file));
G = computeGeometry(G);
layers = [1 2 8 12 13];
K = logNormLayers(G.cartDims, [50 300 100 20]*milli*darcy, ...
                  'sz', [51 3 3], 'indices', layers);
rock = makeRock(G, bsxfun(@times, K*milli*darcy,[1 1 .1]), 0.2);

W = verticalWell([], G, rock, G.cartDims(1), 1, [], 'InnerProduct', 'ip_tpf');
W = verticalWell(W, G, rock, 14, 36, [],'InnerProduct', 'ip_tpf');

clf
plotCellData(G, log10(rock.perm(:,1)),'EdgeAlpha',.1); view(3); 
plotWell(G, W, 'height', 30,'FontSize',12); 
colormap(jet), axis tight off, view(250, 50)


%% Uniform/METIS partition
% We consider two possible partitions. The first is a standard
% load-balanced partition in index space. The second is an unstructured
% partitioning by the METIS graph library that adapts to the underlying
% geology. To this end, we transmissibilities as edge-weights in the
% graph-partitioning algorithm of METIS so that it tries to make grid
% blocks having as homogeneous permeability as possible.
if useMetis
    hT = computeTrans(G, rock);
    p0 = partitionMETIS(G, hT, 7*9*4, 'useLog', true);
else
    p0 = partitionUI(G, [7 9 1]); max(p0)
end
figure, myPlotPartition(G, W, p0, 400);

%% Split fault faces
% It may also be advantageous to split blocks across faults, even if this
% will introduce some very small blocks
if useMetis
    % hT(G.faces.tag==1) = 1e-10*min(hT);
    % pf = partitionMETIS(G, hT, 7*9*4, 'useLog', true, 'ufactor',5);
    pf = p0;
else
    pf = processPartition(G, p0, find(G.faces.tag==1)); max(pf)
end
figure, myPlotPartition(G, W, pf, 400);

%% Radial refinement
% The function 'refineNearWell' takes a set of points and partitions these
% according to the distance in the xy-plane from a single well point. Here,
% we will this function to refine the coarse blocks that contains wells.
% The first well is placed in the corner and we partition the perforated
% well block into five radial sections. The width of the radial sections is
% set to decay as log(r). For the second well, we refine all the
% neighboring blocks surrounding the well block using a number of angular
% sectors that increases as we move radially out from the well point.
angSectors = {1, [1,4,6,9]};
radSectors = [4 4];
pw = pf;
for i=1:numel(W)
    wc    = W(i).cells(1);
    wpt   = G.cells.centroids(wc,:);
    pwv   = pw(wc);
    if i==1
        % Pick all cells inside the well block
        cells = (pw==pwv);
    else
        % Use adjacency matrix to compute nearest neighbors. We start with
        % a vector e which is equal one in the well block and zero
        % elsewhere. Multiplying by the adjacency matrix A will set the
        % value in each block equal the sum over the block and its
        % face-neighbors. After two multiplications, all blocks surrounding
        % the initial block should have value larger than one.
        CG    = generateCoarseGrid(G, pw);
        A     = getConnectivityMatrix(getNeighbourship(CG),true,CG.cells.num);
        rblk  = zeros(CG.cells.num,1); rblk(pwv)=1; rblk((A*A*rblk)>1) = 1;
        cells = rblk(pw)>0;
    end
    pts   = G.cells.centroids(cells,:);
    out   = refineNearWell(pts, wpt, 'angleBins', angSectors{i}, ...
        'radiusBins', radSectors(i), 'logbins', true, 'maxRadius', inf);
    pw(cells) = max(pw) + out;
end
pw = compressPartition(pw);
figure, myPlotPartition(G, W, pw, 400);

%% Refine partition in the vertical direction
% For the Cartesian partition, we end by refining the partition in the
% vertical direction so that it follows the layering in the petrophysical
% data, which we assume is known.
if ~useMetis
    pK = rldecode((1:numel(layers)-1)',diff(layers));
    [~,~,k]=gridLogicalIndices(G);
    pk = compressPartition((pK(k)-1)*max(pw)+pw);
    pk = processPartition(G, pk);
    figure, myPlotPartition(G, W, pk, 400);
end