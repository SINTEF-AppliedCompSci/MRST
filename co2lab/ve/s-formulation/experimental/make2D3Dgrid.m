function [g,  region3D, grdecl] = make2D3Dgrid(grdecl,region3D, varargin)
%cartDims = [10,1,2]
%grdecl = simpleGrdecl(cartDims, physDims, 0, 'undisturbed', true);


opt = struct('g3D', [], 'g_top', []);
opt = merge_options(opt, varargin{:});

g3D = opt.g3D;
g_top = opt.g_top;

if isempty(g3D)
g3D = processGRDECL(grdecl);
g3D = computeGeometry(g3D(1));
end
if isempty(g_top)
   g_top = topSurfaceGrid(g3D);
end
[ijk{1:3}] = ind2sub(g3D.cartDims, g3D.cells.indexMap(:));
ijk        = [ijk{:}];

% make 2D/3D grid
% reshape with cartDims zcoorn(:,:,1) = layer 1, zcoorn(:,:,2:3) = layer 2
zcorn=reshape(grdecl.ZCORN,2*grdecl.cartDims);

% i, j, index to cells that should be changed
ijRegion = ijk(~region3D,:);
ijRegion = sortrows(ijRegion, [2 1 3]);
[ij, ix_first] = unique(ijRegion(:,1:2), 'rows', 'first');
[ij, ix_last] = unique(ijRegion(:,1:2), 'rows', 'last');

i_ix = ij(:,1); %Region2D(:,1); 
j_ix = ij(:,2); %Region2D(:,2); 
k_start = ijRegion(ix_first, 3);
k_stop = ijRegion(ix_last, 3);
   
figure;

plotGrid(g3D, 'faceColor', 'none');
plotGrid(g3D,  find(~region3D));

ACTNUM = zeros(numel(grdecl.ACTNUM),1);


%for all cells that should be 2D/VE   
for c = 1:numel(i_ix)
i_change = unique(mcolon(max(1, 2*i_ix(c)-1), min(2*max(ijk(:,1)), 2*i_ix(c))));
j_change = unique(mcolon(max(1, 2*j_ix(c)-1), min(2*max(ijk(:,2)), 2*j_ix(c))));

if k_start(c) == k_stop(c), continue, end

k_change = unique(mcolon(max(1, 2*k_start(c)-1), 2*k_stop(c)));
%change all layers except top and bottom
%for iz=2:grdecl.cartDims(3)*2; % leave top coordinate, set other coordinate equal to bottom coordinate 
  zcorn(i_change, j_change, k_change(2:end-1)) = ...
     repmat(zcorn(i_change, j_change, k_change(end)), [1, 1, numel(k_change)-2]);
  removed = ijRegion(ix_first+1:ix_last, :); % [repmat([i_ix(c), j_ix(c)], [], k_stop(c)-k_start(c)), (2:grdecl.cartDims(3))'];
  removed = sub2ind(g3D.cartDims, removed(:,1), removed(:,2), removed(:,3));
  grdecl.ACTNUM(removed) = 0;
  %ACTNUM(removed) = 1;
%end
end
%ACTNUM = reshape(ACTNUM, g3D.cartDims);

grdecl.ZCORN=zcorn(:);
g=processGRDECL(grdecl);
g = computeGeometry(g(1));

%find number to 3D cells in original grid
region3D_all = zeros(numel(grdecl.ACTNUM),1);
region3D_all(g3D.cells.indexMap(region3D)) = find(region3D);

region3DCells = region3D_all(g.cells.indexMap);
region3D = region3DCells ~=0;

cells3D_all = zeros(numel(grdecl.ACTNUM),1);
cells3D_all(g3D.cells.indexMap) = 1:g3D.cells.num;

cells3D = cells3D_all(g.cells.indexMap);

clf
plotGrid(g, 'faceColor', 'none');
plotGrid(g, find(region3D));
view(3)


% Find mapping to topsurfacegrid
topcells3D = g_top.columns.cells(g_top.cells.columnPos(1:end-1));
[trash, mapTopSurfaceCell] = ismember(cells3D, topcells3D);

% cell number of cells in topSurfaceGrid
g.cells.mapTopSurface = mapTopSurfaceCell;
% region information
g.cells.region3D = region3D;
% cell numbering of cells in the original 3D grid
g.cells.inx3D = cells3D;


% add information from the topSurface grid3D grid, only for top cells
g.cells.columnPos = g_top.cells.columnPos;
g.columns = g_top.columns;
g.cells.H = g_top.cells.H;
g.cells.z = g_top.cells.z;

% add info about boundary between 2D and 3D grid
isInt = ~any(g.faces.neighbors==0,2);

faces2D = any(~region3D(g.faces.neighbors(isInt,:)),2);
faces3D = any(region3D(g.faces.neighbors(isInt,:)),2);
faceNum = find(isInt);
bndFaces2D3D  = faceNum(faces2D & faces3D);
g.facesBnd.index = bndFaces2D3D;

bndCells = reshape(g.faces.neighbors(bndFaces2D3D,:)',[],1);

%unique because one 2D cell has several 3D neighbors
bndCells2D = unique(bndCells(~region3D(bndCells)));

bndCells3D = unique(bndCells(region3D(bndCells)));

% cell and cellface info for faces on boundary between 2D and 3D
g.facesBnd.cells2D = bndCells2D;
g.facesBnd.cells3D = bndCells3D;
faceIx2D = mcolon(g.cells.facePos(bndCells2D),g.cells.facePos(bndCells2D+1)-1).';
faceIx3D = mcolon(g.cells.facePos(bndCells3D),g.cells.facePos(bndCells3D+1)-1).';
g.facesBnd.cellFace2D = faceIx2D(ismember(g.cells.faces(faceIx2D), bndFaces2D3D));
g.facesBnd.cellFace3D = faceIx3D(ismember(g.cells.faces(faceIx3D), bndFaces2D3D));

cellIx = mcolon(g.cells.columnPos(mapTopSurfaceCell(bndCells2D)), ...
g.cells.columnPos(mapTopSurfaceCell(bndCells2D)+1)-1).';

% for each face, 3D cellindex in underlying 3D grid 
g.facesBnd.map3D = g.columns.cells(cellIx);  
g.facesBnd.centroid3D = g3D.cells.centroids(g.columns.cells(cellIx),:);

% include the topgrid and 3D grid
g.topSurfaceGrid = g_top;
g.parent         = g3D;
end
