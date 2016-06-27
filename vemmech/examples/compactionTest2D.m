%% Define parameters
griddim=3;
opt=struct('L',[1000 1000 100],...
            'cartDims',ones(1,griddim)*ceil((1e3).^(1/griddim)),...
            'grid_type','triangle',...
            'disturb',0.0,... %parameter for disturbing grid
            'E',0.3*1e9,...  %youngs modolo
            'nu',0.4,...
            'islinear',false);% poiso ratio'
opt.L=[15 15 3];% sbed model
opt.hanging=false;
opt.free_top=true;
opt.triangulate=true;% triangulate some faces
opt.vertical=false;%only valid for norne
opt.gravity_load=true;
opt.top_load=false;
opt.gtol=0.1e-1;
opt.ref=10;
opt.pressure=false;
opt.mycase='cartgrid';
opt.twist=true;opt.disturb=0.05;
opt.force_method='cell_force';
%opt.force_method='dual_grad_type';
%opt.mycase='pebi';
opt.cartDims=[10 10];
opt.L=[15*10/10 15]
opt.flipgrid=false;
%
mycase='cartgrid'

G=squareGrid(opt.cartDims,opt.L,'grid_type',opt.mycase);
if(opt.triangulate)
    face=unique(G.cells.faces(any(bsxfun(@eq,G.cells.faces(:,2),[3,4]),2),1));
    face=unique(G.cells.faces(any(bsxfun(@eq,G.cells.faces(:,2),[1,4]),2),1));
    G=triangulateFaces(G,face)
    G=sortEdges(G);
end
if(opt.twist)
G=twister(G,opt.disturb)
end
if(opt.flipgrid)
   G=flipGrid(G);
end
G=mrstGridWithFullMappings(G);
G=computeGeometryCalc(G);
%
figure()
clf,plotGrid(G,'Marker','*')
%%
[el_bc,load] =makeCompactionTest(G,opt);

%% define rock parameters
Ev=repmat(opt.E,G.cells.num,1);nuv=repmat(opt.nu,G.cells.num,1);
C=Enu2C(Ev,nuv,G);
% solve system also with profiler on
profile off;profile on
lsolve=@mldivide;
if(strcmp(mycase,'box'))
lsolve=@agmg;
end
bbsize = 30000-(G.griddim-2)*20000;
G=computeGeometry(G);
%G=computeGeometryCalc(G);
%G.weights.cell_nodes=wc;
profile off;profile on;
%
% change load for pressure since it is a constant vector
minp=-max(sum(load(G.faces.centroids).*G.faces.centroids,2));

if(opt.pressure)
    pressure =@(X) -sum(load(X).*X,2)-minp;
    load =@(X) load(X)*0;
    %{
    bccoord=G.faces.centroids(el_bc.force_bc.faces,:);
    press_boundary=sum(pressure(bccoord));
    bsign=(2* (G.faces.neighbors(el_bc.force_bc.faces,2)==0)-1);
    normals = -bsxfun(@times,G.faces.normals(el_bc.force_bc.faces,:),bsign);
    el_bc.force_bc.force=el_bc.force_bc.force+bsxfun(@times,normals,press_boundary);
    %}
else
   pressure =@(X) -sum(load(X).*X,2)*0.0;  
end
%%
uu=VEM_linElast(G,C,el_bc,load,'linsolve',lsolve,'blocksize',bbsize,'force_method',opt.force_method,'pressure',pressure(G.cells.centroids));
profile off;profile viewer

%%
div=VEM_div(G);
fac=1;
%clf,plotNodeDataDeformed(G,uu(:,1),fac*uu,'EdgeColor','none');colorbar
%clf,plotCellDataDeformed(G,div*reshape(uu',[],1),fac*uu,'EdgeColor','none');colorbar
figure(),clf
%plotNodeDataDeformed(G,uu(:,3),fac*uu,'EdgeColor','none');colorbar
plotNodeDataDeformed(G,uu(:,G.griddim),fac*uu);cb=colorbar;%view(3)
cdiv=div*reshape(uu',[],1);
%plotCellDataDeformed(G,cdiv,fac*uu,'EdgeColor','none');colorbar
%plotFaces(G,bc{1}.face,'FaceColor','r')
%plotFaces(G,bc{2}.face,'FaceColor','b')
%set(cb,'YTick',[-2e-2:1e-2:0])
%axis off
%title(mycase)
%%
if(~isempty(el_bc.force_bc))
    ff=abs(el_bc.force_bc.force(1,G.griddim));
else
    ff=0;
end
start=max(G.faces.centroids(:,G.griddim));
top=min(G.faces.centroids(:,G.griddim));

%ana=@(z) ff*(z-start)./(lambda+mu);
%ana=@(z) ff*(z-start)./(C(1,1))+double(opt.gravity_load)*10*300*((z).^2-(start).^2)/C(1,1);
%if(G.griddim==2)
%fac=1*300;
%else
    fac=100*300/2;
%end
ana=@(z) ff*(z-start)./(C(1,1))-double(opt.gravity_load)*fac*((top-start).^2 - (z-top).^2)/C(1,1);
diva=@(z) (ff./C(1,1))-double(opt.gravity_load)*fac*(-2*(z-top))/C(1,1);
%
z=G.nodes.coords(:,G.griddim);
z(abs(ana(z))<max(abs(ana(z)))*1e-2)=nan;
zl=unique(z);
figure(),
subplot(4,1,1)
plot(z,uu(:,G.griddim),'*',zl,ana(zl))
subplot(4,1,2)
err=(uu(:,G.griddim)-ana(z))./max(abs(ana(z)));
plot(z,err,'*')
subplot(4,1,3)
div=VEM_div(G);
plot(G.cells.centroids(:,G.griddim),div*reshape(uu',[],1)./G.cells.volumes,'*');
subplot(4,1,4)
div=VEM_div(G);
zc=G.cells.centroids(:,G.griddim);
diverr=(div*reshape(uu',[],1)./G.cells.volumes-diva(zc))./max(abs(diva(zc)));
plot(zc,diverr,'*');

%%
figure(),
subplot(2,1,1) 
%err=diverr;
nodes=find(abs(err)>1e-3);
icells=mcolon(G.nodes.cellPos(nodes),G.nodes.cellPos(nodes+1)-1);
cells=unique(G.nodes.cells(icells));
[~,nodes]=max(abs(err));
icells=mcolon(G.nodes.cellPos(nodes),G.nodes.cellPos(nodes+1)-1);
mcells=unique(G.nodes.cells(icells));
plotGrid(G,cells)
plotGrid(G,mcells,'FaceColor','r')
subplot(2,1,2)
cells=find(abs(diverr)>1e-3);
plotGrid(G,cells)
[~,mcells]=max(abs(diverr));
plotGrid(G,mcells,'FaceColor','r')

plotGrid(G,cells)
