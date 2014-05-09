%% Construct Top-Surface Grid of the Johansen Formation
% The Johansen formation is a candidate site for large-scale CO2 storage
% offshore the south-west coast of Norway. We consider a sector model that
% has been constructed based on available seismic and well data and stored
% in the Eclipse input format (GRDECL). The model has five vertical layers
% in the Johansen formation and five shale layers above and one below in
% the Dunhil and Amundsen formations. The shale layers are removed and we
% construct a 2D grid of the top surface, assuming that the major
% fault is sealing. Moreover, we identify all outer boundaries that are
% open to flow. A more thorough presentation of the geological model can be
% found in the script <matlab:edit('showJohansen.m') showJohansen.m>
%
% The grid and rock structured constructed in the following can be used for
% subsequent 3D and/or VE simulations of CO2 injection and migration and
% are therefore stored to file to avoid time-consuming processing.
%
% The data files necessary to run the example can be downloaded from the 
% <http://www.sintef.no/Projectweb/MatMorA/Downloads/Johansen/ MatMoRA
% website>.

%% Load model and construct VE grid
% Load grid geometry - you will most likely have to change the path,
% depending upon where you have stored the Johansen data-set
clear, clc
sector = fullfile(ROOTDIR, 'data', 'johansen', 'NPD5');
fprintf('    Reading %s\n', sector);
grdecl = readGRDECL([sector, '.grdecl']);

% Load permeability and porosity
K = reshape(load([sector, '_Permeability.txt'])', [], 1);
p = reshape(load([sector, '_Porosity.txt'])',     [], 1);
grdecl.PERMX=K;grdecl.PERMY=K;grdecl.PERMZ=0.1*K;
grdecl.PORO=p;
clear sector;

% Remove shale layers Dunhil and Amundsen
disp('    Constructing grids ...');
grdecl.ACTNUM(K.*milli*darcy<0.11 *milli*darcy) = 0;

% Construct grid structure. By removing Dunhil and Amundsen, the grid will
% consist of multiple components. We choose the largest one as our
% reservoir.
G = processGRDECL(grdecl, 'verbose', false);
G = computeGeometry(G(1));
%clear grdecl

% Construct VE grid
[G_top, G] = topSurfaceGrid(G);

grdecl.ACTNUM(:) = 0;
grdecl.ACTNUM(G.cells.indexMap) = 1;

[ijk{1:3}] = ind2sub(G.cartDims, G.cells.indexMap(:));
ijk        = [ijk{:}];

% etter plume spread
region3D = (ijk(:,1) > 43 & ijk(:,1) < 53  );
region3D = region3D & (ijk(:,2)> 46 & ijk(:,2) < 56  );

%region3D = false(G.cells.num,1);

% figure;
% plotGrid(G, 'faceColor', 'none');
% plotGrid(G, find(region3D));

[G_c, trash, region3D] = make2D3Dgrid(grdecl,region3D, 'g3D', G, 'g_top', G_top);

clear grdecl
% Construct structure with petrophyiscal data. Vertical permeability is set
% to 0.1 times the horizontal permeability. NB!
rock.perm = bsxfun(@times, [1 1 0.1], K(G.cells.indexMap)).*milli*darcy;
rock.poro = p(G.cells.indexMap);


%%{
rock.perm(:) = median(rock.perm(:,1));
rock.poro(:) = median(rock.poro);
% nB: should use 3D perm in 3D area and 2D perm in 2D area.. 
rock_c.perm = rock.perm(1)*ones(G_c.cells.num,1);
rock_c.poro = rock.poro(1)*ones(G_c.cells.num,1);

rock2D    = averageRock(rock, G_top);


clear p K;

%% FIND FAULT - set inner boundary
% The main fault is assumed to be sealing and must therefore be represented
% as an inner boundary in the 2D grid.  Start by locating cells on each
% side of fault, which is defined to be between cells with index i=43 and
% i=44 for j<44.
cells2D_1 = find(G_top.cells.ij(:,1) == 43 & G_top.cells.ij(:,2) <= 44);
cells2D_2 = find(G_top.cells.ij(:,1) == 44 & G_top.cells.ij(:,2) <= 44);

% Plot the cells on opposite sides
figure;
plotGrid(G_top, 'faceColor', 'none');
plotGrid(G_top, cells2D_1, 'faceColor', 'r')
plotGrid(G_top, cells2D_2, 'faceColor', 'g')
axis tight off,
title('Cells on opposite sides of sealing fault'),

% Find the faces at the fault. Construct a mapping facesMat defined such
% that facesMat(i,j)=k if face <k> is shared by cells <i> and <j>.
facesMat = sparse(double(G_top.faces.neighbors(:,1))+1, ...
   double(G_top.faces.neighbors(:,2))+1, 1:G_top.faces.num);
facesMat = facesMat + facesMat';
faultFaces2D = diag( facesMat(cells2D_1+1, cells2D_2+1) );

% Make internal boundary and compute the geometry of the resulting grid
G_top = makeInternalBoundary(G_top, faultFaces2D);
G_top = computeGeometryVE(G_top);

% The function 'topSurfaceGrid' does not handle faults correctly when the
% cells on opposite sides are not phyiscally in contact with each other.
% Instead of producing a discontinuity in the 'z' value, an average value
% is used.  Hence, we need to manually reset the 'z' value of these cells
% (marked in red and green in the plot) to avoid an incorrect flow 
G_top.cells.z([cells2D_1; cells2D_2]) = ...
   G.faces.centroids(G_top.cells.map3DFace([cells2D_1; cells2D_2]), 3);
clear cells2D_1 cells2D_2 facesMat faultFaces2D


%% FIND PRESSURE BOUNDARY
% Setting boundary conditions is unfortunately a manual process and may
% require some fiddling with indices, as shown in the code below. Here, we
% identify the part of the outer boundary that is open, i.e., not in
% contact with one of the shales (Dunhil or Amundsen).

% boundary 3D
nx = G.cartDims(1); ny=G.cartDims(2); nz=G.cartDims(3);
ix1 = boundaryFaceIndices(G, 'BACK', 1:nx-6, 1:4, 1:nz);
ix2 = boundaryFaceIndices(G, 'LEFT', 1:20,   1:ny, 1:nz);
ix3 = boundaryFaceIndices(G, 'RIGHT', 1:nx, ny-10:ny, 1:nz);
ix4 = boundaryFaceIndices(G, 'FRONT', 1:nx/2-8, ny/2:ny, 1:nz);

figure;
subplot(1,2,1)
plotGrid(G, 'faceColor', 'none', 'EdgeAlpha', 0.1)
plotFaces(G, ix1, 'r');
plotFaces(G, ix2, 'g');
plotFaces(G, ix3, 'y')
plotFaces(G, ix4, 'm');
view(-7,60), axis tight off;
title('3D grid');

bcIx = [ix1; ix2; ix3; ix4];
clear ix1 ix2 ix3 ix4


% boundary 3D
nx = G_c.cartDims(1); ny=G_c.cartDims(2); nz=G_c.cartDims(3);
ix1 = boundaryFaceIndices(G_c, 'BACK', 1:nx-6, 1:4, 1:nz);
ix2 = boundaryFaceIndices(G_c, 'LEFT', 1:20,   1:ny, 1:nz);
ix3 = boundaryFaceIndices(G_c, 'RIGHT', 1:nx, ny-10:ny, 1:nz);
ix4 = boundaryFaceIndices(G_c, 'FRONT', 1:nx/2-8, ny/2:ny, 1:nz);

figure;
subplot(1,2,1)
plotGrid(G_c, 'faceColor', 'none', 'EdgeAlpha', 0.1)
plotFaces(G_c, ix1, 'r');
plotFaces(G_c, ix2, 'g');
plotFaces(G_c, ix3, 'y')
plotFaces(G_c, ix4, 'm');
view(-7,60), axis tight off;
title('3D/2D grid');

bcIx_c = [ix1; ix2; ix3; ix4];
clear ix1 ix2 ix3 ix4

% boundary 2D
nx = G_top.cartDims(1); ny=G_top.cartDims(2);
ix1 = boundaryFaceIndices(G_top, 'BACK',  1:nx-6, 1:4, []);
ix2 = boundaryFaceIndices(G_top, 'LEFT',  1:20, 1:ny,  []);
ix3 = boundaryFaceIndices(G_top, 'RIGHT', 1:nx, ny-10:ny, []);
ix4 = boundaryFaceIndices(G_top, 'FRONT', 1:nx/2-8, ny/2:ny, []);

%remove faces connected to main fault
ix1 = ix1(G_top.faces.centroids(ix1,2)>6.714*1e6);
ix2 = ix2(G_top.faces.centroids(ix2,1)>5.4*1e5);
%
subplot(1,2,2)
plotGrid(G_top, 'faceColor', 'none', 'EdgeAlpha', 0.1)
plotGrid(G_top, sum(G_top.faces.neighbors(ix1,:),2), 'faceColor', 'r')
plotGrid(G_top, sum(G_top.faces.neighbors(ix2,:),2), 'faceColor', 'g')
plotGrid(G_top, sum(G_top.faces.neighbors(ix3,:),2), 'faceColor', 'y')
plotGrid(G_top, sum(G_top.faces.neighbors(ix4,:),2), 'faceColor', 'm')
axis tight off
title('2D grid of top surface');

bcIxVE = [ix1; ix2; ix3; ix4];
clear ix1 ix2 ix3 ix4 nx ny nz

%% Store data
disp('    Writing JohansenVEgridSIAM.mat')
save JohansenVEgridSIAM.mat G G_top G_c region3D rock rock_c rock2D bcIx bcIx_c bcIxVE
