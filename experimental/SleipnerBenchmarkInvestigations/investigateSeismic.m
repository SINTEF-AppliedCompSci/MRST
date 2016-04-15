%addpath([pwd,'/helpers/'])
%mrstModule add co2lab-plain
files={'data/2008-CO2Anomaly-mincurv.asc',...
    'data/1994-TopLayer9-kriging.asc',...
    'data/2008-StructuralThickness_SandLayer9.asc'}
%files={'data/1994-TopLayer9-kriging.asc',...
%    'data/1994-TopLayer9-kriging_depth.asc'}
%files={'data/1994-TopLayer9-kriging.asc',...
%    'data/2008-TopLayer9-kriging.asc'}
%files={'data/2008-CO2Anomaly-mincurv.asc',...
%    'data/2008-TopLayer9-kriging.asc'}

for k=1:numel(files)
%%
%a = importdata('data/1994-topLayer9-kriging.asc',' ',20)
%a = importdata('data/1994-topLayer9-kriging_Depth.asc',' ',20)
%a = importdata('data/2008-BottomCO2Anomaly-mincurv_Depth.asc',' ',20)
%a = importdata('data/2008-CO2Anomaly-mincurv.asc',' ',20)
a = importdata(files{k},' ',20);
% cartdim linje
%%
lin=a.textdata{14}
b=strsplit(lin,':')
assert(strcmp(b{1},'# Grid_size'));
tmp=strsplit(b{end},'x')
assert(numel(tmp)==2);
grid_size=nan(2,1);
for i=1:2
    grid_size(i)=str2num(tmp{i});
end
lin=a.textdata{15}
b=strsplit(lin,':')
assert(strcmp(b{1},'# Grid_space'))
tmp=strsplit(b{end},',')
assert(numel(tmp)==4);
grid_space=nan(4,1);
for i=1:numel(tmp)
    grid_space(i)=str2num(tmp{i});
end
val=cell(3,1);
for i=1:numel(val)
    %%
   val{i} =nan(grid_size');
   ind=sub2ind(size(val{i}),a.data(:,4),a.data(:,5));
   val{i}(ind)=a.data(:,i);
end
datas{k}.val=val;
end

%%
%figure(),hold on
for k=1:numel(datas)
    figure()
    mesh(datas{k}.val{1},datas{k}.val{2},datas{k}.val{3})
end
%%
%figure()
%for k=1:numel(datas)
%    mesh(datas{k}.val{1},datas{k}.val{2},datas{k}.val{3})
%end
%%
k1=2;k2=1;
figure()
val=nan(size(datas{k2}.val{1}));
ind=isfinite(datas{k2}.val{1});
val(ind)=interp2(datas{k1}.val{1}',datas{k1}.val{2}',datas{k1}.val{3}',datas{k2}.val{1}(ind),datas{k2}.val{2}(ind));
Ft=@(x,y) interp2(datas{k1}.val{1}',datas{k1}.val{2}',datas{k1}.val{3}',x,y);
%F=makegriddedinterp(datas{k1}.val{1},datas{k1}.val{2},datas{k1}.val{3});
%val(ind)=F(datas{k2}.val{1}(ind),datas{k2}.val{2}(ind));
dval=datas{k2}.val{3}-val;
mesh(datas{k2}.val{1},datas{k2}.val{2},-dval);
%%

figure(33),clf
subplot(1,2,1)
title([files{k1},' - ',files{k2}])
%mesh(datas{k2}.val{1},datas{k2}.val{2},-dval);view([-1 -2 5])
pcolor(datas{k2}.val{1},datas{k2}.val{2},-dval),shading flat,caxis([0 15]),colorbar
subplot(1,2,2)
k=3;
title(files{k})
%mesh(datas{k}.val{1},datas{k}.val{2},datas{k}.val{3});view([-1 -2 5])
pcolor(datas{k}.val{1},datas{k}.val{2},datas{k}.val{3}),shading flat,caxis([0 15]),colorbar
%%
for i=1:numel(datas)
   F{i}=@(x,y) interp2(datas{k}.val{1}',datas{k}.val{2}',datas{k}.val{3}',x,y)
end

return
%%
[ G, Gt, rock, rock2D ] = makeSleipnerModelGrid()
%%
%
%deckfile='../../co2lab/data/sleipner/SLEIPNER_DECK.DATA';
%grdecl=readGRDECL(deckfile);
%%
%G1=processGRDECL(grdecl);
%%
% [G,Gt]  = makeSleipnerVEmodel()
figure(),clf
k=2
mesh(datas{k}.val{1},datas{k}.val{2},datas{k}.val{3})
plotGrid(G,'EdgeAlpha',0.1)

%%
figure()
%plotCellData(Gt,Ft(Gt.cells.centroids(:,1),Gt.cells.centroids(:,2))),colorbar
plotCellData(Gt,Ft(Gt.cells.centroids(:,1),Gt.cells.centroids(:,2))-Gt.cells.z),colorbar

%%
tdata= readSeismicAsc('data/1994-TopLayer9-kriging_depth.asc');
%Hdata= readSeismicAsc('data/2008-StructuralThickness_SandLayer9.asc');
%%
X=reshape(Gt.cells.centroids(:,1),Gt.cartDims);
Y=reshape(Gt.cells.centroids(:,2),Gt.cartDims);
H=reshape(Gt.cells.H,Gt.cartDims);
%{
x=x(:,end:-1:1);
y=y(:,end:-1:1);
H=H(:,end:-1:1);
FH=@(x,y) interp2(x,y,H,x,y);
%}
FH=scatteredInterpolant(X(:),Y(:),H(:));
Fttmp=scatteredInterpolant(Gt.cells.centroids(:,1),Gt.cells.centroids(:,2),Gt.cells.z);
Ft=@(vec) Fttmp(vec(:,1),vec(:,2));
%Gn = ascData2Grid(tdata,'cartdims',[100,100],'FH',FH)
Gn = ascData2Grid(tdata,'cartdims',Gt.cartDims,'FH',FH)
Gtn=topSurfaceGrid(computeGeometry(Gn))
Ftntmp=scatteredInterpolant(Gtn.cells.centroids(:,1),Gtn.cells.centroids(:,2),Gtn.cells.z);
Ftn=@(vec) Ftntmp(vec(:,1),vec(:,2));
%%
% compare co2 plume
figure(),
plotGrid(Gt,'FaceColor','r')
zt=Ft(Gtn.cells.centroids);
ztn=Ftn(Gtn.cells.centroids);
plotCellData(Gtn,ztn-zt),colorbar
%5

