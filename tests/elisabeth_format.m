
%%
%%{
clear all
tic
[simp_grid,meta_simp,petro] = getAtlasGrid('Utsirafm', 'coarsening', 1,'refining',7,'make_deck',false);%,'refining',refining);
clear meta_simp;
petro=petro{1}
simp_grid=simp_grid{1};
cell_thick=(simp_grid.thick(1:end-1,1:end-1)+simp_grid.thick(2:end,1:end-1)+simp_grid.thick(2:end,2:end)+simp_grid.thick(1:end-1,2:end))/4;
cell_top=(simp_grid.top(1:end-1,1:end-1)+simp_grid.top(2:end,1:end-1)+simp_grid.top(2:end,2:end)+simp_grid.top(1:end-1,2:end))/4;
%
cell_X=(simp_grid.X(1:end-1,1:end-1)+simp_grid.X(2:end,1:end-1)+simp_grid.X(2:end,2:end)+simp_grid.X(1:end-1,2:end))/4;
cell_Y=(simp_grid.Y(1:end-1,1:end-1)+simp_grid.Y(2:end,1:end-1)+simp_grid.Y(2:end,2:end)+simp_grid.Y(1:end-1,2:end))/4;
%
centroid=1.0e+06*[0.4386, 6.4750];%-simp_grid.orig;
no_cell=isnan(cell_thick) | isnan(cell_top);
cell_thick(no_cell)=-100;
cell_top(no_cell)=-100;
cartDims=simp_grid.nodeDims-1;
prod(cartDims)
toc;
h_init=10*exp(- ((cell_X-centroid(1)).^2 +  (cell_Y-centroid(2)).^2)/1e4.^2);
pcolor(cell_X,cell_Y, 100*h_init+cell_top ) ;shading flat;colorbar;

%%
tic;
datadir = '/home/hnil/tmp/';

fileId = fopen(fullfile(datadir,'top.bin'), 'w');

fprintf(fileId, '%0.7f\n', cell_top(:));
fclose(fileId);
%
fileId = fopen(fullfile(datadir,'resHeight.bin'), 'w');
fprintf(fileId, '%0.7f\n', cell_thick(:));
fclose(fileId);
%
dx =  simp_grid.dx;%(Gt.nodes.coords(2) - Gt.nodes.coords(1));
dim = [cartDims(1), cartDims(2), dx];
fileId = fopen(fullfile(datadir,'dimensions.txt'), 'w');
fprintf(fileId, '%d\n', dim);
fclose(fileId);
clear dim;

fileId = fopen(fullfile(datadir,'co2height.bin'), 'w');
fprintf(fileId, '%0.7f\n', h_init(:));
toc;


%%
return
%%
mesh(simp_grid.X,simp_grid.Y,simp_grid.thick)
set(gca,'ZDir','reverse')
%%

require deckformat
%%
clear all
profile off
profile on

refining=1;
coarsening = 10;
[utsira_grdecl,utsira_meta] = getAtlasGrid('Utsirafm', 'coarsening', coarsening);%,'refining',refining);
utsira_grdecl = utsira_grdecl{1};
%
% Process and make top surface grid
tic;
utsira_grdecl=refineGrdecl(utsira_grdecl,[1 1 1])
G3D_utsira = processGRDECL(utsira_grdecl);
clear utsira_grdecl;
Gtop_utsira = topSurfaceGrid(G3D_utsira);
toc
%profile off
%profile viewer
Gtop_utsira.cells.num/1e6

flux_on_grid=ones(Gtop_utsira.faces.num,1);
flux_all_cart=zeros(prod(Gtop_utsira.cartDims),1);
flux_all_cart(Gtop_utsira.faces.cartFaceMap)=flux_on_grid;
profile off
profile viewer



info_top=utsira_meta_simp{2};
info_top=refineInfo(info_top,1);
Gt = processAAIGrid(info_top.meta, info_top.data, true, [1, 1]);
plotCellData(Gt,Gt.cells.z);
toc
%}

%%
datasets = readAtlasGrids('Utsirafm', 10);
