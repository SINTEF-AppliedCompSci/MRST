function [g, region3D] = addGtopFields(g, g_top, g3D, region3D)

%find number to 3D cells in original grid
region3D_all = zeros(prod(g3D.cartDims),1);
region3D_all(g3D.cells.indexMap(region3D)) = find(region3D);

region3DCells = region3D_all(g.cells.indexMap);
region3D = region3DCells ~=0;

cells3D_all = zeros(prod(g3D.cartDims),1);
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
