%% See how the sleipner field is placed on the utsira formation
%  When considering CO2 storage the Sleipner case is the most famous
%  example. Here we see how model desribed in 
%  is placed on the Utsira formation. We first plot the Utsira formation
% get utsira grdecl and meta data

require deckformat

coarsening = 5;
[utsira_grdecl,utsira_meta] = getAtlasGrid('Utsirafm', 'coarsening', coarsening);
utsira_grdecl = utsira_grdecl{1};

% Process and make top surface grid
G3D_utsira = processGRDECL(utsira_grdecl);
Gtop_utsira = topSurfaceGrid(G3D_utsira);

%% Get sleipner data
deck_sleipner=readEclipseDeck(fullfile(VEROOTDIR,'data', 'sleipner', 'SLEIPNER_DECK.DATA'));
% 
%grdecl=deck.GRID;
%grdecl.ACTNUM(grdecl.PERMX<200) = 0;
% this model layers 3:36 is not shale layer
grdecl_sleipner=deck_sleipner.GRID;
%G3D_sleipner_fine=processGRDECL(grdecl);
% cut out the shale part of the model
grdecl_cc=cutGrdecl(grdecl_sleipner,...
    [1,grdecl_sleipner.cartDims(1);1,grdecl_sleipner.cartDims(2);3,36]);
% make a one layer model
grdecl_cc=coarseGrdecl(grdecl_cc,[1 1,grdecl_cc.cartDims(3)]);
coord_tmp=reshape(grdecl_cc.COORD,3,[])';
coord=coord_tmp;
% do mapaxes explicitely
coord(:,1:2)=mapAxes(coord_tmp(:,1:2), deck_sleipner.GRID.MAPAXES);
coord=coord';
grdecl_mapax=grdecl_cc;grdecl_mapax.COORD=coord(:);
clear grdecl_cc grdecl_tmp
G3D_sleipner = processGRDECL(grdecl_mapax);
G3D_sleipner = computeGeometry(G3D_sleipner);
Gtop_sleipner=topSurfaceGrid(G3D_sleipner);

%% plot 3D models
figure(1),clf,
plotGrid(G3D_utsira,'FaceColor','none','edgea', .05, 'edgec', 'b');
plotGrid(G3D_sleipner,'FaceColor','r','EdgeColor','none');
Gtop_utsira_surf=processAAIGrid(utsira_meta{2}.meta, utsira_meta{2}.data, true, [2 2]);
plotGrid(Gtop_utsira_surf,'FaceColor','none','edgea', .05, 'edgec', 'k');
axis off;axis tight;view([-1 -1 3])
%%
figure(2),clf
plotCellData(Gtop_utsira, Gtop_utsira.cells.z, 'facea', .6, 'edgea', .05, 'edgec', 'k');
%plotGrid(Gtop_utsira,'FaceColor','none')
set(gcf,'Color',[0.8 0.8 0.8]);set(gca,'Color',[0.8 0.8 0.8]);
set(gca,'LineWidth',1);colorbar;axis equal;axis tight;
plotGrid(Gtop_sleipner, 'FaceColor','r','EdgeColor','none')
%set(gcf,'Color',[0 0 0])
view(2)
box on
drawnow;

%% Plot just a small portion of utsira
% find region around sleipner with a fine utsira model
figure(3),clf
center=sum(Gtop_sleipner.cells.centroids,1)/Gtop_sleipner.cells.num;
[d,cind_s] = min(sum(bsxfun(@minus,Gtop_sleipner.cells.centroids,center).^2,2));%#ok
z_sleipner_center=Gtop_sleipner.cells.z(cind_s);
coord_min=min(Gtop_sleipner.nodes.coords);
coord_max=max(Gtop_sleipner.nodes.coords);
ndiag=norm(coord_min-coord_max);
dist_center=sum(bsxfun(@minus,Gtop_utsira.cells.centroids(:,1:2),center).^2,2);
p=sqrt(dist_center)<4*ndiag;
%Gtop_utsira_small=processPartition(Gtop_utsira,p);
plotCellData(Gtop_utsira, Gtop_utsira.cells.z, find(p), 'facea', .6, 'edgea', .05, 'edgec', 'k');
set(gcf,'Color',[0.8 0.8 0.8]);set(gca,'Color',[0.8 0.8 0.8]);
set(gca,'LineWidth',1);colorbar;axis tight;
plotCellData(Gtop_sleipner, Gtop_sleipner.cells.z, 'facea', .6, 'edgea', .05, 'edgec', 'k');
%set(gcf,'Color',[0 0 0])
axis on; view([-1 -1 3])


%% move grid to match
figure(4),clf
[d,cind_utsira]=min(dist_center);%#ok
z_utsira_center=Gtop_utsira.cells.z(cind_utsira);

p=sqrt(dist_center)<2*ndiag;
%Gtop_utsira_small=processPartition(Gtop_utsira,p);
plotCellData(Gtop_utsira, Gtop_utsira.cells.z, find(p), 'facea', 0.6, 'edgea', .05, 'edgec', 'k');
set(gcf,'Color',[0.8 0.8 0.8]);set(gca,'Color',[0.8 0.8 0.8]);
set(gca,'LineWidth',1);colorbar;axis tight;
Gtop_sleipner_mod=Gtop_sleipner;
z_shift=z_utsira_center-z_sleipner_center;
Gtop_sleipner_mod.nodes.z=Gtop_sleipner_mod.nodes.z+z_shift;
Gtop_sleipner_mod.cells.z=Gtop_sleipner_mod.cells.z+z_shift;
disp(['Shift surface with ', num2str(z_shift),' meter']); 
plotCellData(Gtop_sleipner_mod, Gtop_sleipner_mod.cells.z, 'facea', .8, 'edgea', .05, 'edgec', 'k');
%set(gcf,'Color',[0 0 0])
axis on; view([-1 -1 3])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find gravity streamlines on utsira from sleipner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
require streamlines


dist_center=sum(bsxfun(@minus,Gtop_utsira.cells.centroids(:,1:2),center).^2,2);

[d,cind_utsira]=min(dist_center);
internal=all(Gtop_utsira.faces.neighbors>0,2);
flux=zeros(Gtop_utsira.faces.num,1);
flux(internal)=Gtop_utsira.cells.z(Gtop_utsira.faces.neighbors(internal,2))-Gtop_utsira.cells.z(Gtop_utsira.faces.neighbors(internal,1));
state=struct('flux',flux);
%
nc=ceil(Gtop_utsira.cells.num/500);
s_cells=[1:nc:Gtop_utsira.cells.num]';%#ok
[S,T,C]=pollock(Gtop_utsira,state,[s_cells, repmat([0.5,0.5],numel(s_cells),1)],'reverse',true);
%[S,T,C]=pollock(Gtop_utsira,state);%,Gtop_utsira.cells.centroids(cind_utsira,:))
%S = reshape([S, repmat({[nan, nan]}, [numel(S),1])]',[], 1);
%
y=linspace(0,utsira_meta{2}.meta.cellsize*(utsira_meta{2}.meta.ncols-1),utsira_meta{2}.meta.ncols)...
    +utsira_meta{2}.meta.yllcorner;
x=linspace(0,utsira_meta{2}.meta.cellsize*(utsira_meta{2}.meta.nrows-1),utsira_meta{2}.meta.nrows)...
    +utsira_meta{2}.meta.xllcorner;
[X,Y]=meshgrid(x,y);
data=abs(utsira_meta{2}.data);
data=data';
%data=data(end:-1:1,:);
data(data==0)=nan;
% make coarse representation for conveninience
nc=1;X=X(1:nc:end,1:nc:end);Y=Y(1:nc:end,1:nc:end);data=data(1:nc:end,1:nc:end);
Z=cell(size(S));
for k=1:numel(S)
   Z{k}=interp2(X,Y,data, S{k}(:,1), S{k}(:,2));
end


%%
figure(5),clf
% use the original meta data to surf make the surface
%surf(X,Y,data+10),shading interp
surf(X,Y,data+10),shading interp
% plot the part inside a small ratidus of the sleipner data
%plotGrid(Gtop_utsira,'FaceColor','r')
plotGrid(Gtop_utsira,sqrt(dist_center)<0.27*ndiag,'FaceColor','r')
%%
hold on
for k=1:numel(S)
    plot3(S{k}(:,1), S{k}(:,2), Z{k}-3,'k-','LineWidth',1);
end
hold off

%plotGrid(Gtop_sleipner,'FaceColor','r')
set(gca,'Zdir','reverse')
% plot streamlines for upward moving plume for some cells in a small radius
s_cells=find(sqrt(dist_center)<2*ndiag);
s_cells=s_cells(1:ceil(numel(s_cells)/20):numel(s_cells));
[SN,TN,CN]=pollock(Gtop_utsira,state, [s_cells, repmat([0.5,0.5],numel(s_cells),1)],...
    'reverse',true,'maxsteps', 5000);%#ok
ZN=cell(size(SN));
% move the streamline to 3D
for k=1:numel(SN)
   ZN{k}=interp2(X,Y,data, SN{k}(:,1), SN{k}(:,2));
end
hold on;
% plot streamlines
for k=1:numel(SN)
    plot3(SN{k}(:,1), SN{k}(:,2), ZN{k}-1,'r-','LineWidth',2);
end
% make streamlines downwards on th e utsira formation for from cell cind_utsira
% starting at the midt point
[SN,TN,CN]=pollock(Gtop_utsira,state, [cind_utsira, 0.5,0.5],'reverse',false,'maxsteps', 5000);
ZN=cell(size(SN));
% move the streamline to 3D
for k=1:numel(SN)
   ZN{k}=interp2(X,Y,data, SN{k}(:,1), SN{k}(:,2));
end
hold on;
% plot streamlines (acturally only one)
for k=1:numel(SN)
    plot3(SN{k}(:,1), SN{k}(:,2), ZN{k}-5,'r-','LineWidth',2);
end
axis off;axis tight
set(gca,'Zdir','reverse')
view([1 -4 6])
% find trapping structure
trap_str=findTrappingStructure(Gtop_utsira);
Gtop_utsira_flat=trap_str.Gtop;
z_spill_loc=trap_str.z_spill_loc;

%
plotGrid(Gtop_utsira,z_spill_loc>0,'EdgeColor','none')
drawnow;


