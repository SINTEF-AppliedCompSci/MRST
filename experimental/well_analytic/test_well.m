%figure(1),clf
%
H=10;
mycases={};
mycases{1}=struct('L',[11000,11000,H],'dim',[[1 1]*10+1,1],'well_case','vertw');
mycases{2}=struct('L',[11000,11000,H],'dim',[[1 1]*100+1,10],'well_case','single_h');
mycases{3}=struct('L',[11000,11000,H],'dim',[[1 1]*40+1,10],'well_case','double_h');
mycases{4}=struct('L',[11000,11000,H],'dim',[[1 1]*100+1,1],'well_case','vertw');
mycases{5}=struct('L',[11000,11000,H],'dim',[[1 1]*100+1,1],'well_case','single_h');
%mycases{5}=struct('L',[11000,11000,H],'dim',[[1 1]*40*3+1,10],'well_case','point');

%L=[11000,11000,100]
%dim=[[1 1]*40+1,10];
%well_case='double_h';

%{
grdecl=simpleGrdecl(dim,0,'physDim',L,'flat',true,'undisturbed',true);

grdecl=refineGrdecl(grdecl
G=processGRDECL(grdecl);
%}
sols={};
%%
for k=1:numel(mycases)
mycase=mycases{k};
well_case=mycase.well_case;
L=mycases{k}.L;
dim=mycases{k}.dim;

Lw=1000;
dx=L(1)/(2*dim(1));
G=cartGrid(dim,L);
G=computeGeometry(G);

switch well_case
    case 'single_h'
        distx=sqrt(sum((G.cells.centroids(:,1)-L(1)/2).^2,2));
        [md,i] = min(distx);
        [I,J,K]=ind2sub(G.cartDims,G.cells.indexMap);
        cells=find((distx<dx) & abs(G.cells.centroids(:,2)-L(2)/2)<Lw/2 & K==G.cartDims(3));
        wdir='x'
    case 'double_h'
        %%
        figure()
        wdis=Lw/2
        distx1=sqrt(sum((G.cells.centroids(:,1)-((L(1)/2+wdis))).^2,2));
        distx2=sqrt(sum((G.cells.centroids(:,1)-((L(1)/2-wdis))).^2,2));
        distx=min(distx1,distx2);
        [md,i] = min(distx);
        [I,J,K]=ind2sub(G.cartDims,G.cells.indexMap);
        cells=find((distx<dx) & abs(G.cells.centroids(:,2)-L(2)/2)<Lw/2 & K==G.cartDims(3));
        wdir='x'
    case 'point'
        distxy=sqrt(sum(bsxfun(@minus,G.cells.centroids(:,1:2),L(1:2)/2).^2,2));
        [md,i] = min(distxy);
        [I,J,K]=ind2sub(G.cartDims,G.cells.indexMap);
        cells=find((distxy<dx) & K==G.cartDims(3));
        wdir='z';
    case 'vertw'
        distxy=sqrt(sum(bsxfun(@minus,G.cells.centroids(:,1:2),L(1:2)/2).^2,2));
        [md,i] = min(distxy);
        %[I,J,K]=ind2sub(G.cartDims,G.cells.indexMap);
        cells=find((distxy<dx) );
        wdir='z';
        
        
        
end
clf
  plotGrid(G,'FaceColor','none')
  plotGrid(G,cells,'FaceColor','r')       
%%
rock=struct('perm',darcy*ones(G.cells.num,1),'poro',ones(G.cells.num,1));
gravity off
fluid=initSingleFluid('mu',1*centi*poise,'rho',1000);
W = addWell([], G, rock, cells, 'Radius',0.3,'Dir',wdir,'Val',1e6/year,'Type','rate');

sides={'Xmin','Xmax','Ymin','Ymax'}
faces=[];
for i=1:numel(sides)
    faces=[faces;boundaryFaceIndices(G,sides{i})];
end
bc=addBC([],faces,'pressure',0*barsa);

T=computeTrans(G,rock);
state0=initState(G,W,100*barsa);
state = incompTPFA(state0, G, T, fluid, 'bc', bc, 'src', [], ...
                        'wells', W);
%%
X3D=reshape(G.cells.centroids(:,1),G.cartDims);
Y3D=reshape(G.cells.centroids(:,2),G.cartDims)
if(G.cartDims(3)>1)
X=squeeze(X3D(:,:,1));
Y=squeeze(Y3D(:,:,2));
else
X=X3D;
Y=Y3D; 
end
p3D=reshape(state.pressure,G.cartDims);
figure()
mesh(X,Y,squeeze(p3D(:,:,1))/barsa)
%sols{k}=struct('X',X,'Y',Y,'p3D',p3D,'cells',cells,'G',G,'state',state)
sols{k}=struct('cells',cells,'G',G,'state',state,'W',W)
end

%%
%figure()
%mesh(X,Y,squeeze(p3D(:,:,1))/barsa)
%hold on
%mesh(X,Y,(squeeze(p3D(:,:,1))-squeeze(p3D(:,:,end)))/barsa)
% lines at top perpendicular to well
%cc=colormap();
col={'r','b','g','m','k'}
leg={};
figure(),clf,hold on
for k=1:numel(sols)
    sol=sols{k};
    leg{k}=mycases{k}.well_case;
    state=sols{k}.state;
    G=sol.G;
    
    X3D=reshape(G.cells.centroids(:,1),G.cartDims);
    Y3D=reshape(G.cells.centroids(:,2),G.cartDims);
    if(G.cartDims(3)>1)
        X=squeeze(X3D(:,:,1));
        Y=squeeze(Y3D(:,:,2));
    else
        X=X3D;
        Y=Y3D;
    end
    p3D=reshape(state.pressure,G.cartDims);
    if(true)
        nx=ceil(G.cartDims(1)/2);
    assert(all(X(nx,:)==L(1)/2))
    yc=Y(nx,:);
    pp=squeeze(p3D(nx,:,1));
    plot(yc,pp/barsa,col{k})
    else
            ny=ceil(G.cartDims(2)/2);
    assert(all(Y(:,ny)==L(2)/2))
    yc=X(:,ny);
    pp=squeeze(p3D(:,ny,1));
    plot(yc,pp/barsa,col{k})
    end
end
ax=axis();
re=sols{1}.W.re;
%line([re,re]+L(1)/2,[ax(3) ax(4)])
%line(-[re,re]+L(1)/2,[ax(3) ax(4)])
%line([Lw/2,Lw/2]+L(1)/2,[ax(3) ax(4)])
%line(-[Lw/2,Lw/2]+L(1)/2,[ax(3) ax(4)])

legend(leg{:})
%%
rad=5500;
q=1e6/year;
%logl =@(Z,lw) -real(2*log(1/lw*(-(sq1(Z-lw)+sq2(Z+lw)))));
logl =@(Z,lw) -real(2*log(1/lw*((sq1(Z-lw)+sq2(Z+lw)))));
fana=@(x) -((q*centi*poise)/(2*pi*1*darcy*H)).*log(x/rad);
fanav=@(x) -((q*centi*poise)/(2*pi*1*darcy*H)).*logl(x/rad,Lw/(2*rad))
plot(yc,fana(yc-rad)/barsa,yc,fanav(yc-rad)/barsa)



