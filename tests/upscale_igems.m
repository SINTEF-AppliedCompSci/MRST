% to do trap analysis we need  the following modules
mrstModule add internal/vertical-equil
mrstModule add coarsegrid
mrstModule('add','/home/hnil/heim/SVN/simmatlab/projects/statoil_upscaling')
mrstModule('add','/home/hnil/heim/SVN/simmatlab/branches/halvor/mrst_extra/gridprocessing')
%addpath('/home/hnil/heim/SVN/simmatlab/branches/halvor/mrst_extra/gridprocessing')
mycase = 'igems_grdecl'
mycase = 'sara_3d'
%mycase = 'igems_surf'
for hh=2:2
switch mycase
    case 'sara_3d'
        %%
        if(hh==1)
          a=load('/home/hnil/heim/SVN/simmatlab/branches/mrst-reorg/modules/internal/vertical-equil/examples/igems/data/sara/g_top_FMM.mat');
          b=load('/data/igems/30times60km/matlab_airbus/igems_original_FMM_2zones_z1_uniform.mat');          
          z_bottum=reshape(b.Gt.cells.z+b.Gt.cells.H,b.Gt.cartDims)';          
          Gt=a.g_top_FMM;
          Gt.cells.H=z_bottum(:)-Gt.cells.z(:);
          Gt.cells.H=Gt.cells.H-mean(Gt.cells.H)+100;
        else
          a=load('/home/hnil/heim/SVN/simmatlab/branches/mrst-reorg/modules/internal/vertical-equil/examples/igems/data/sara/g_top_OSS.mat');
          b=load('/data/igems/30times60km/matlab_airbus/igems_original_OSS_2zones_z1_uniform.mat');
          z_bottum=reshape(b.Gt.cells.z+b.Gt.cells.H,b.Gt.cartDims);          
          %Gt=a.g_top_OSS;
          Gt=b.Gt;
          Gt.cells.z=a.g_top_OSS.cells.z;
          Gt.faces.z=a.g_top_OSS.faces.z;
          Gt.cells.H=z_bottum(:)-Gt.cells.z(:);
          Gt.cells.H=Gt.cells.H-mean(Gt.cells.H)+100;          
        end
        clear b
        clear a
        Gt.griddim=2;
        Gt.cells.cellNodes=getCellNodes(Gt);

        
        %G=a.G;
        %Gt=a.Gt;
        %Gt.cells.z=-Gt.cells.z;
        %Gt.faces.z=-Gt.faces.z;
        clear a;
        gravity([0 0 9.8])
        %newcartDims=Gt.cartDims(1:2)
        if(false)
            %newcartDims=[200 200];
            %[I,J]=meshgrid(50:50+newcartDims(1),200:200+newcartDims(1));
            newcartDims=[100 100];
            [I,J]=meshgrid(100:100+newcartDims(1),200:200+newcartDims(1));
            I=I';J=J';
            newcartDims=newcartDims+1;
            cells=sub2ind(Gt.cartDims,I(:),J(:));
        else
            nn=4
            newcartDims=  ([25 25]*nn);
            p = partitionCartGrid(Gt.cartDims, floor([Gt.cartDims(1:2)./newcartDims]));
            cells=find(p==11);
        end
        Gt_old=Gt;
        %Gt=extractSubgrid(Gt,find(p==1));
        %[Gt, gcells, gfaces, gnodes] = extractSubgrid(Gt, find(p==1))
        %[Gt, gcells, gfaces, gnodes] = extractSubgrid(Gt, find(p==11));        
        [Gt, gcells, gfaces, gnodes] = extractSubgrid(Gt_old,cells);
        Gt=computeGeometry(Gt);
        Gt.cells.z=Gt_old.cells.z(gcells);
        Gt.cells.H=Gt_old.cells.H(gcells);
        %assumme bottom periodic
        
        
        %Gt.cells.H=100*ones(numel(gcells),1);
        
        Gt.faces.z=Gt_old.faces.z(gfaces);
        Gt.nodes.z=Gt_old.nodes.z(gnodes);
        Gt.cells.indexMap=[1:Gt.cells.num];
        Gt.cartDims=newcartDims;
        Gt.cells.cellNodes=getCellNodes(Gt);
        rock.poro=0.25*ones(Gt.cells.num,1);
        rock.perm=ones(Gt.cells.num,1).*Gt.cells.H*darcy;
        
    case 'igems_grdecl'
        if(hh==1)
        a=load('/data/igems/30times60km/matlab_airbus/igems_original_FMM_2zones_z1_uniform.mat')
        else
        a=load('/data/igems/30times60km/matlab_airbus/igems_original_OSS_2zones_z1_uniform.mat')
        end
        G=a.G;
        Gt=a.Gt;
        %Gt.cells.z=-Gt.cells.z;
        %Gt.faces.z=-Gt.faces.z;
        clear a;
        nn=4
        newcartDims=  ([25 25]*nn);
        gravity([0 0 9.8])
        %newcartDims=Gt.cartDims(1:2)
        p = partitionCartGrid(Gt.cartDims, [Gt.cartDims(1:2)./newcartDims]);
        Gt_old=Gt;
        %Gt=extractSubgrid(Gt,find(p==1));
        %[Gt, gcells, gfaces, gnodes] = extractSubgrid(Gt, find(p==1))
        [Gt, gcells, gfaces, gnodes] = extractSubgrid(Gt, find(p==11));
        Gt=computeGeometry(Gt);
        Gt.cells.z=Gt_old.cells.z(gcells);
        Gt.cells.H=Gt_old.cells.H(gcells);
        %assumme bottom periodic
        
        
        %Gt.cells.H=100*ones(numel(gcells),1);
        
        Gt.faces.z=Gt_old.faces.z(gfaces);
        Gt.nodes.z=Gt_old.nodes.z(gnodes);
        Gt.cells.indexMap=[1:Gt.cells.num];
        Gt.cartDims=newcartDims;
        Gt.cells.cellNodes=getCellNodes(Gt);
        rock.poro=0.25*ones(Gt.cells.num,1);
        rock.perm=ones(Gt.cells.num,1).*Gt.cells.H*darcy;
    case 'igems_surf'
        gravity([0 0 9.8])
        % we can load surfaces from the igems data set the possible names are
        % the directories  in data/surfaces
        name='flatNP2';
        name='OSS';
        % there is 100 different surfaces we choose number
        i=2;
        % For testing it is usefull to use a coarse variant
        coarse=[4,4];
        %[Gt,rock2D,SVE,rock,G]=readIGEMSIRAP(name,i,coarse)
        [Gt,rock2D,rock,G]=readIGEMSIRAP(name,i,coarse);
        if(false)
            
            p = partitionCartGrid(G.cartDims, [Gt.cartDims(1:2)/25,1]);
            G=extractSubgrid(G,find(p==1));
            %G.cartDims=
            G= computeGeometry(G);
            Gt=topSurfaceGrid(G);
            Gt=computeGeometry(Gt);
        else
            nn=4
            newcartDims=  ([25 25]*nn)./coarse;
            %newcartDims=Gt.cartDims(1:2)
            p = partitionCartGrid(Gt.cartDims, [Gt.cartDims(1:2)./newcartDims]);
            Gt_old=Gt;
            %Gt=extractSubgrid(Gt,find(p==1));
            %[Gt, gcells, gfaces, gnodes] = extractSubgrid(Gt, find(p==1))
            [Gt, gcells, gfaces, gnodes] = extractSubgrid(Gt_old, find(p==(4/nn)^2*10))
            Gt=computeGeometry(Gt);
            Gt.cells.z=Gt_old.cells.z(gcells);
            %Gt.cells.H=Gt_old.cells.H(gcells);
            Gt.cells.H=100*ones(numel(gcells),1);
            
            Gt.faces.z=Gt_old.faces.z(gfaces);
            Gt.nodes.z=Gt_old.nodes.z(gnodes);
            Gt.cells.indexMap=[1:Gt.cells.num];
            Gt.cartDims=newcartDims;
        end
        %%
        g_vec=[0         0    9.8066];
        figure(1),clf
        %plotCellData(Gt,p);
        %
        
        
        
        rock.poro=ones(Gt.cells.num,1);
        rock.perm=ones(Gt.cells.num,1).*Gt.cells.H*darcy;
% define periodic boundary
        %return
    case 'test_case'
        g_vec=[0         0    9.8066];
        % simple test grid for finding traps
        L=100;L_p=L/2;H=10;
        G=cartGrid([10 10 1],[L L H]);
        ff=@(coord) sin(pi*coord(:,1)/L_p)*1+sin(pi*coord(:,1)/L_p).*(sign(sin(pi*coord(:,1)/L_p))+1).*sin(pi*coord(:,2)/L_p)*0.1+coord(:,1)*0.01;
        G.nodes.coords(:,3)=ff(G.nodes.coords)+G.nodes.coords(:,3);
        %100+G.nodes.coords(:,3)+G.nodes.coords(:,1)*0.0...
        %    -4*sin(pi*G.nodes.coords(:,1)/L_p);
        G=computeGeometry(G);
        Gt=topSurfaceGrid(G);
        % overide 2d surface calulalations
        Gt=computeGeometry(Gt);
        Gt.cells.z=ff(Gt.cells.centroids);
        Gt.faces.z=ff(Gt.faces.centroids);
        Gt.cells.H=ones(Gt.cells.num,1)*100;
        Gt.grav_pressure=@(g,omega) zeros(numel(g.cells.faces(:,1)),1);
        rock.poro=ones(Gt.cells.num,1);
        rock.perm=ones(Gt.cells.num,1).*Gt.cells.H*darcy;
% define periodic boundary
    otherwise
        error('No such case')
end
%%
figure(44)
subplot(2,1,1)
cla,plotCellData(Gt_old,Gt_old.cells.z);colorbar
axis equal
%plotGrid(Gt_old,find(p==11),'FaceColor','none','EdgeAlpha',0.1)
plotGrid(Gt_old,cells,'FaceColor','none','EdgeAlpha',0.1)
subplot(2,1,2)
cla,mesh(reshape(Gt.cells.centroids(:,1),Gt.cartDims),reshape(Gt.cells.centroids(:,2),Gt.cartDims),reshape(Gt.cells.z,Gt.cartDims))
hold on
mesh(reshape(Gt.cells.centroids(:,1),Gt.cartDims),reshape(Gt.cells.centroids(:,2),Gt.cartDims),reshape(Gt.cells.z+Gt.cells.H,Gt.cartDims))
figure(33),clf,mesh(reshape(Gt.cells.z,Gt.cartDims))
%return
%cla,plotGrid(Gt)
%%
%% Start upscaling process

clear G;
G=Gt;
clear Gt;

%%


%%
%return
if(false)
  cellno=rldecode([1:G.cells.num]',diff(G.cells.facePos));
dhfz=-(G.cells.z(cellno)-G.faces.z(G.cells.faces(:,1)));  

N = G.faces.neighbors;
i = ~any(N==0, 2);
dzf=zeros(G.faces.num,1);
dzf(i)=G.cells.z(N(i,2)) - G.cells.z(N(i,1));
sgn=2*(G.faces.neighbors(G.cells.faces(:,1),1)==cellno)-1;
dzf_tmp=accumarray(G.cells.faces(:,1),dhfz.*sgn);
[dzf,dzf_tmp]
end
%return
%% plot traps
%%


dp_scale=1;
if(G.griddim>1)
   bcl{1}=pside([],G,'LEFT',0);dp{1}=0;
   bcr{1}=pside([],G,'RIGHT',0);
   bcl{2}=pside([],G,'BACK',0);dp{2}=0;
   bcr{2}=pside([],G,'FRONT',0);
else
   bcl{1}=struct('face',1,'value',dp_scale);dp{1}=dp_scale;
   bcr{1}=struct('face',G.faces.num,'value',dp_scale);
end
if(G.griddim>2)
   bcl{3}=pside([],G,'TOP',0);dp{3}=0;
   bcr{3}=pside([],G,'BOTTOM',0);
end
% find grid dimensions of the upscaling domain
L=nan(1,G.griddim);
for i=1:G.griddim
   L(i) = max(G.faces.centroids(:,i))-min(G.faces.centroids(:,i));
end
% for all the fluids are define outside for experimenting
fluid_pure=initSingleFluid('mu',1,'rho',1);
%fluid_org = initSimpleVEFluidSForm('mu',[0.05 0.5]*centi*poise,'rho',[700 1020],'sr',[0 0],'height',G.cells.H);
fluid_org = initSimpleVEFluidSForm('mu',[6e-5 6.9e-4],'rho',[700 1020],'sr',[0 0],'height',G.cells.H);
fluid_nc=fluid_org;
fluid_nc.pc =@(state) zeros(size(state.s,1),1);
% set input
fluids{1}=fluid_pure;fluids{2}=fluid_org;fluids{3}=fluid_nc;
% define upscale values
dSS=1/(25*5);
dS=1/10;
sat_vec=[dS:dS:1-4*dS];%1-dS];
sat_vec=[dS:dS:0.2];%1-dS];
sat_vec=[dS:dS:1-4*dS];%1-dS];
%sat_vec=[0.5];%1-dS];
sat_vec=[dSS:dSS:0.3,0.3+dS:dS:1-dS]
dp_vec=1;initcap=false;
[sat_mat,kr,perm,Kkr,creport]=relperm_upscale_periodic_gravity(G,rock,fluids,dp_scale,sat_vec,[],[], initcap);
figure(1)
%%
figure(1),clf
plot(G.cells.centroids(:,1),G.cells.z,'b*',...    
     G.cells.centroids(:,1),G.cells.z-G.cells.H,'b*')      
     disp('hei')
figure(2)
plot(sat_mat,Kkr{1}(:,:),'*-',sat_mat,Kkr{2}(:,:),'*-')
plot(sat_mat,Kkr{1}(:,:),'*-')%,sat_mat,Kkr{2}(:,1),'*-')%,sat_mat,Kkr{2}(:,1),'*-')
%%
kr={};
for i=1:2
kr{i}=bsxfun(@rdivide,Kkr{i}(:,:),[perm(:)]');
kr{i}=kr{i}(:,[1,4]);
end
figure(3)
plot(sat_mat,kr{1},sat_mat,kr{2})
if(hh==1)
  save(['relperm_all_fmm'],'sat_mat','kr')  
else
  save(['relperm_all_oss'],'sat_mat','kr')  
end
save(['relperm_all',num2str(hh)],'sat_mat','kr')
end

return
%%
figure(1)
clf,
mesh(reshape(Gt_old.cells.z,Gt_old.cartDims)')
hold on; 
mesh(reshape(Gt_old.cells.z+Gt_old.cells.H,Gt_old.cartDims)');
hold on
mesh(reshape(a.Gt.cells.z+a.Gt.cells.H,a.Gt.cartDims));
hold on
mesh(reshape(a.Gt.cells.z,a.Gt.cartDims));
mesh(reshape(G.cells.z,G.cartDims)',10*ones(G.cartDims))
hold on; 
mesh(reshape(G.cells.z+G.cells.H,G.cartDims)',10*ones(G.cartDims));