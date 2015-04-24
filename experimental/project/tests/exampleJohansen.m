%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See how different data sets of johansen is related
%
nc=3;
[johansenfm_grdecl,johansenfm_meta]=getAtlasGrid('Johansenfm','coarsening',nc)
johansenfm_grdecl=johansenfm_grdecl{1};
johansenfm_meta=johansenfm_meta{1};
% get johansen small model
jdir = fullfile(VEROOTDIR, 'data', 'johansen');
sector = 'NPD5';
fprintf('    Reading %s\n', sector);
sector = fullfile(jdir, sector);
johansen_grdecl = readGRDECL([sector '.grdecl']);
johansen_large_grdecl=readGRDECL(...
    fullfile(VEROOTDIR, 'data', 'johansen','FULLFIELD_IMAXJMAX.GRDECL')...
    );
K = reshape(load([sector, '_Permeability.txt'])', [], 1);
p = reshape(load([sector, '_Porosity.txt'])',     [], 1);
johansen_grdecl.PERMX=K;johansen_grdecl.PERMY=K;johansen_grdecl.PERMZ=0.1*K;
johansen_grdecl.PORO=p;
%%
if(true)
 johansen_grdecl=cutGrdecl(johansen_grdecl,[10 50; 11 40; 1 11]);
    %johansen_grdecl= coarseGrdecl(johansen_grdecl,[nc nc,1]);  
%% get johansen large model
%johansen_large_grdecl=cutGrdecl(johansen_large_grdecl,[70 110; 50 100; 10 14]);
johansen_large_grdecl=cutGrdecl(johansen_large_grdecl,[40 140; 40 140; 1 johansen_large_grdecl.cartDims(3)]);
%johansen_large_grdecl=cutGrdecl(johansen_large_grdecl,[1 149; 1 189; 10 14]);
end

G3D_johansenfm=processGRDECL(johansenfm_grdecl);
G3D_johansen=processGRDECL(johansen_grdecl,'SplitDisconnected',false);
G3D_johansen_large=processGRDECL(johansen_large_grdecl,'SplitDisconnected',false);
figure(1),clf
plotGrid(G3D_johansenfm,'FaceColor','b','edgea',0.2,'facea',0.2);
plotGrid(G3D_johansen,'FaceColor','r','edgea',0.2,'facea',0.2);
plotGrid(G3D_johansen_large,'FaceColor','none','edgea',0.2,'facea',0.2);
axis off;axis tight;
%{
plotGrid(G3D_johansenfm,'FaceColor','b','edgea',0.2,'facea',0.2);
plotCellData(G3D_johansen,johansen_grdecl.PERMX(johansen_grdecl.ACTNUM>0))
plotCellData(G3D_johansen_large,johansen_large_grdecl.PERMX(johansen_large_grdecl.ACTNUM>0))
%}
colorbar
%nc=5;
%johansen_grdecl=cutGrdecl(johansen_grdecl,[10 50; 11 40; 1 11]);
%johansen_grdecl= coarseGrdecl(johansen_grdecl,[nc nc,1]);
%% get johansen large model
%johansen_large_grdecl=cutGrdecl(johansen_large_grdecl,[70 110; 50 100; 10 14]);
%johansen_large_grdecl=cutGrdecl(johansen_large_grdecl,[40 140; 40 140; 10 14]);
%hist(log10(johansen_large_grdecl.PERMX(johansen_large_grdecl.ACTNUM>0 & johansen_large_grdecl.PERMX>0)))
%johansen_large_grdecl= coarseGrdecl(johansen_large_grdecl,[nc nc,1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The data set the surface given as the johansen formation in from NPD in
% the co2 atlas is the surface on top of the johansen data when permeability 
% less than  0.11 milli*darcy is removed. For the johansen_large data it is
% the top of the lowest layer which is separated from the rest by a gap in
% the data.
%johansen_grdecl.ACTNUM(K.*milli*darcy<0.11 *milli*darcy) = 0;
johansen_grdecl.ACTNUM(johansen_grdecl.PERMX<0.11) = 0;
johansen_large_grdecl=cutGrdecl(johansen_large_grdecl,[1 johansen_large_grdecl.cartDims(1); 1 johansen_large_grdecl.cartDims(2); 10 14]);
johansen_large_grdecl.ACTNUM(johansen_large_grdecl.PERMX<100)=0;


G3D_johansen=processGRDECL(johansen_grdecl,'SplitDisconnected',false);
G3D_johansen_large=processGRDECL(johansen_large_grdecl,'SplitDisconnected',false);
%% plot matching 3D formations
figure(2),clf,
plotGrid(G3D_johansenfm,'FaceColor','b','edgea',0.2,'facea',0.2);
plotGrid(G3D_johansen,'FaceColor','r','edgea',0.2,'facea',0.2);
plotGrid(G3D_johansen_large,'FaceColor','none','edgea',0.2,'facea',0.2);
axis off;axis tight;
% generate topsurface grids
Gtop_johansenfm=topSurfaceGrid(G3D_johansenfm);
Gtop_johansen_large=topSurfaceGrid(G3D_johansen_large);
Gtop_johansen=topSurfaceGrid(G3D_johansen);
% plot matching surfaces
figure(3),clf
plotGrid(Gtop_johansenfm,'FaceColor','b','edgea',0.2,'facea',0.2);
plotGrid(Gtop_johansen,'FaceColor','r','edgea',0.2,'facea',0.2);
plotGrid(Gtop_johansen_large,'FaceColor','none','edgea',0.2,'facea',0.2);
axis off;axis tight;