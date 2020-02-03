mrstModule add mimetic incomp

%Convergence tests
clear all;
close all;
%clc

%mrstVerbose on
mrstVerbose off

%% 2D test: rotating anisotropic permeability, u=sin(pi*x)*sin(pi*y)+1
Lx=1;Ly=1;
kesi=1e-3; % anisotropy ratio mild 1e-1, strong 1e-3
GG=4;
%n=[8 16 32]'% 64]'% 128 256]';
n=2.^[3:7];
nc=zeros(numel(n),1);
ep=zeros(numel(n),2);
ef=zeros(numel(n),2);
picard=zeros(numel(n),2);
hap=zeros(numel(n),1);

tols = [1e-7, 1e-12];

uexact=@(xy) sin(pi*xy(1))*sin(pi*xy(2))+1; %scalar input
uxexact=@(xy) pi*cos(pi*xy(1,:)).*sin(pi*xy(2,:));
uyexact=@(xy) pi*sin(pi*xy(1,:)).*cos(pi*xy(2,:));
fluid=initSingleFluid('mu',1,'rho',1);

for t = 1:numel(tols)
    tol = tols(t);
    for i_con=1:length(n)
        nx=n(i_con);ny=n(i_con);
        switch GG
          case 1
            G=Grid_Normal(Lx,Ly,nx,ny,0.6,'moveYes');
          case 2 % as expected 2 1 conv
            myfile=['./grid_data/tri',num2str(i_con),'.dat'];
            G=importGrid(myfile,'tri');
            G=computeGeometry(G);
          case 3 
            % Mild kesi
            % w/o twister: pressure 2. flux mfd 2, ntpfa 1.5. Almost exact
            % rates for different tols
            % with twister: mfd 2 for both pressure and
            % flux. ntpfa approaches 2 for both pressure and flux
            % Twister doesn't give any centroids outside
            % Cart grid
            G = cartGrid([nx,ny], [Lx,Ly]);
            G = twister(G, 0.1);
            G = computeGeometry(G);
            % Strong kesi: same. pressure and flux -> 2. Lots of
            % picard iterations to reach 1e-12. # iterations
            % basically double with mesh size, starting from 50 => {50,100,200,400,800}
            
          case 4
            % Mild kesi
            % Triangle grid
            % Only first grid has centroids outside convex hull
            % Both pressure and flux approaches 2 for both mfd and ntpfa.
            % No difference for low and high tol.
            G = cartGrid([nx,ny], [Lx,Ly]);
            G = triangleGrid(G.nodes.coords);
            G = twister(G, 0.1);
            G = computeGeometry(G);
            % Strong kesi
            % All grids have fraction outside
            % MFD and ntpfa perform similarly bad.
            % MFD increases up to rate 1.7 for pressure
            % ntpfa increases up to rate 1.4 for pressure
            % MFD increases up to rate 1.12 for flux
            % ntpfa increases up to rate 1.3 for flux
            % Maybe indication of reduced conv rate for ntpfa
        end
        %figure,plotGrid(G,'facecolor','none');axis image;
        nc(i_con)=G.cells.num; clear nx ny

        x=G.cells.centroids(:,1);y=G.cells.centroids(:,2);
        rock.perm=bsxfun(@rdivide,[kesi.*x.^2+y.^2 (kesi-1).*x.*y x.^2+kesi.*y.^2],x.^2+y.^2);
        clear x y
        bc.face=boundaryFaces(G);
        bc.type=repmat({'pressure'},[numel(bc.face),1]);
        bc.value=repmat({uexact},[numel(bc.face),1]);

        x=G.faces.centroids(:,1);y=G.faces.centroids(:,2);
        faceperm=bsxfun(@rdivide,[kesi.*x.^2+y.^2 (kesi-1).*x.*y x.^2+kesi.*y.^2],x.^2+y.^2);
        kxx=faceperm(:,1);kxy=faceperm(:,2);kyy=faceperm(:,3);

        ux=uxexact([x';y'])';
        uy=uyexact([x';y'])';
        vx=-(kxx.*ux+kxy.*uy);vy=-(kxy.*ux+kyy.*uy);
        fa=dot(G.faces.normals,[vx vy],2);
        hf_flux=fa(G.cells.faces(:,1));
        hf2cn=gridCellNo(G);
        sgn=2*(G.faces.neighbors(G.cells.faces(:,1),1)==hf2cn)-1;
        hf_flux=sgn.*hf_flux;
        q=accumarray(hf2cn,hf_flux);
        src=addSource([],1:G.cells.num,q);
        clear x y vx vy hf_flux hf2cn sgn ux uy faceperm kxx kxy kyy

        disp(['Mesh refinement level ' num2str(i_con) ':'])
        %s1=OnePhaseIncompMPFA(G,TransMPFAO(G,rock,bc),src);
        
        % MFD
        xy_faces=G.faces.centroids(bc.face,:);
        %uexact_faces=uexact(xy_faces');
        uexact_faces = zeros(size(xy_faces,1),1);
        for i = 1:size(xy_faces,1)
            uexact_faces(i) = uexact(xy_faces(i,:));
        end
        %uexact_faces=arrayfun(uexact, xy_faces);
        bc_std=addBC([],boundaryFaces(G),'pressure',uexact_faces);
        state0=initResSol(G,1);
        s1=incompMimetic(state0,G,computeMimeticIP(G,rock),fluid,'bc',bc_std,'src',src);

        % NTPFA
        %mrstVerbose on
        [interpFace]=findHAP(G,rock,bc);
        hap(i_con)=interpFace.fraction;
        disp(['fraction of faces with centroids outside convex hull: ', num2str(interpFace.fraction)]);
        interpFace=correctHAP(G,interpFace);
        OSflux=findOSflux(G,rock,bc,interpFace);
        u0=ones(G.cells.num,1);
        s2=FlowNTPFA(G,bc,fluid,[],OSflux,u0,tol,1000,'src',src);
        %s3=PicardNTPFA(G,bc,src,OSflux,u0,tol,500);
        %s3=PicardNMPFA(G,bc,src,OSflux,u0,1e-7,500);
        %picard(i_con,:)=[s2.iter, s3.iter];
        picard(i_con,1)=s2.iter;
        clear u0 OSflux interpFace src bc rock q
        %mrstVerbose off
        
        %x=G.cells.centroids(:,1);y=G.cells.centroids(:,2);
        %ua=uexact(G.cells.centroids');
        ua=zeros(size(G.cells.centroids,1),1);
        for i = 1:size(G.cells.centroids,1)
            ua(i) = uexact(G.cells.centroids(i,:));
        end
        %plot3(x,y,abs(ua-s2.pressure),'.') %,pause
        clear x y
        [ep(i_con,1),ef(i_con,1)]=L2_Norm(G,ua,fa,s1.pressure,s1.flux);
        [ep(i_con,2),ef(i_con,2)]=L2_Norm(G,ua,fa,s2.pressure,s2.flux);
        %[ep(i_con,3),ef(i_con,3)]=L2_Norm(G,ua,fa,s3.pressure,s3.flux);
    end
    disp(tol)
    %plotConvergence(ep,ef,nc,{'MPFA-O','NTPFA','NMPFA'});
    titles = {'MFD',['NTPFA ', num2str(tol)]}
    plotConvergence(ep,ef,nc,titles);
    picard
    %clear fa GG i_con Lx Ly myfile n 
    titles
    disp('Convergence order of pressure solutions:')
    disp(convOrder(nc,ep,2))
    disp('Convergence order of flux solutions:')
    disp(convOrder(nc,ef,2))
    drawnow
end

% %% 3D test with mild anistropy: u(x,y,z)=1+sin(pi*x)*sin(pi*(y+0.5))*sin(pi*(z+1/3))
% lx=1;ly=1;lz=1;
% kxx=1;kxy=0.5;kxz=0;kyy=1;kyz=0.5;kzz=1;
% GG=2;
% n=[5 8 12 18 27]';
% nc=zeros(size(n));
% ep=zeros(length(n),3);
% ef=zeros(length(n),3);
% picard=zeros(numel(n),2);
% hap=zeros(numel(n),1);

% for i_con=1:numel(n)
%     nx=n(i_con);ny=nx;nz=nx;
%     switch GG
%         case 1  % random hexahedras---------------------------------------
%             G=cartGrid([nx ny nz],[lx ly lz]);
%             rn=-0.5+rand(G.nodes.num,3);
%             rn=bsxfun(@times,rn,[lx/nx ly/ny lz/nz]);
%             G.nodes.coords=G.nodes.coords+0.3*rn;
%             G=computeGeometry(G);clear rn
%             figure, plotGrid(G);view(3);axis image;set(gca,'zdir','normal');
%         case 2 % unstructured tetrahedral grid----------------------------
%             myfile=['./grid_data/tet',num2str(i_con),'.dat'];
%             G=importGrid(myfile,'tetBench');
%             G=computeGeometry(G);
%             nc(i_con)=G.cells.num;
%             figure, plotGrid(G);view(3);axis image;set(gca,'zdir','normal');
%     end
%     nc(i_con)=G.cells.num;clear nx ny nz
%     rock.perm=repmat([kxx kxy kxz kyy kyz kzz],G.cells.num,1);
%     bc.face=boundaryFaces(G);
%     bc.type=repmat({'pressure'},[numel(bc.face),1]);
%     bc.value=repmat({@(x)1+sin(pi*x(1))*sin(pi*(x(2)+0.5))*sin(pi*(x(3)+1/3))},[numel(bc.face),1]);
    
%     perm=[kxx kxy kxz;kxy kyy kyz;kxz kyz kzz];
%     x=G.faces.centroids(:,1);y=G.faces.centroids(:,2);z=G.faces.centroids(:,3);
%     grad_u=pi*[cos(pi*x).*sin(pi*(y+0.5)).*sin(pi*(z+1/3)) sin(pi*x).*cos(pi*(y+0.5)).*sin(pi*(z+1/3)) sin(pi*x).*sin(pi*(y+0.5)).*cos(pi*(z+1/3))]';
    
%     velocity=-perm*grad_u;
%     fa=dot(G.faces.normals,velocity',2);
%     hf_flux=fa(G.cells.faces(:,1));
%     hf2cn=gridCellNo(G);
%     sgn=2*(G.faces.neighbors(G.cells.faces(:,1),1)==hf2cn)-1;
%     hf_flux=sgn.*hf_flux;
%     q=accumarray(hf2cn,hf_flux);
%     src=addSource([],1:G.cells.num,q);
%     clear perm x y z grad_u velocity hf_flux sgn q hf2cn
    
%     disp(['Mesh refinement level ' num2str(i_con) ':'])
%     s1=OnePhaseIncompMPFA(G,TransMPFAO(G,rock,bc),src);
    
%     interpFace=findHAP(G,rock,bc);
%     hap(i_con)=interpFace.fraction;
    
%     interpFace=correctHAP(G,interpFace);
%     OSflux=findOSflux(G,rock,bc,interpFace);
%     u0=ones(G.cells.num,1);
%     s2=PicardNTPFA(G,bc,src,OSflux,u0,1e-7,500);
%     s3=PicardNMPFA(G,bc,src,OSflux,u0,1e-7,500);
%     picard(i_con,:)=[s2.iter s3.iter];
%     clear u0 OSflux interpFace src bc rock
    
%     x=G.cells.centroids(:,1);y=G.cells.centroids(:,2);z=G.cells.centroids(:,3);
%     ua=1+sin(pi*x).*sin(pi*(y+0.5)).*sin(pi*(z+1/3));clear x y z
%     [ep(i_con,1),ef(i_con,1)]=L2_Norm(G,ua,fa,s1.pressure,s1.flux);
%     [ep(i_con,2),ef(i_con,2)]=L2_Norm(G,ua,fa,s2.pressure,s2.flux);
%     [ep(i_con,3),ef(i_con,3)]=L2_Norm(G,ua,fa,s3.pressure,s3.flux);
% end
% plotConvergence(ep,ef,nc,{'MPFA-O','NTPFA','NMPFA'});
% clear fa GG i_con Lx Ly myfile n 
% disp('Convergenc order of pressure solutions:')
% disp(convOrder(nc,ep,3))
% disp('Convergenc order of flux solutions:')
% disp(convOrder(nc,ef,3))

% %% 3D test with strong anisotropy: u(x,y,z)=1+sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)
% lx=1;ly=1;lz=1;
% kxx=1;kxy=0;kxz=0;kyy=1;kyz=0;kzz=1e3;
% GG=2;
% n=[5 8 12 18 27]';
% nc=zeros(size(n));
% ep=zeros(length(n),3);
% ef=zeros(length(n),3);
% picard=zeros(numel(n),2);
% hap=zeros(numel(n),1);

% for i_con=1:numel(n)
%     nx=n(i_con);ny=nx;nz=nx;
%     switch GG
%         case 1  % random hexahedras---------------------------------------
%             G=cartGrid([nx ny nz],[lx ly lz]);
%             rn=-0.5+rand(G.nodes.num,3);
%             rn=bsxfun(@times,rn,[lx/nx ly/ny lz/nz]);
%             G.nodes.coords=G.nodes.coords+0.3*rn;
%             G=computeGeometry(G);
%             figure, plotGrid(G);view(3);
%             axis image;set(gca,'zdir','normal');
%             clear x y z ind rn
%         case 2 % unstructured tetrahedral grid----------------------------
%             myfile=['.\grid_data\tet',num2str(i_con),'.dat'];
%             G=importGrid(myfile,'tetBench');
%             G=computeGeometry(G);
%             nc(i_con)=G.cells.num;
%             figure, plotGrid(G);view(3);axis image;set(gca,'zdir','normal');
%     end
%     nc(i_con)=G.cells.num;clear nx ny nz
%     rock.perm=repmat([kxx kxy kxz kyy kyz kzz],G.cells.num,1);

%     bc.face=boundaryFaces(G);
%     bc.type=repmat({'pressure'},[numel(bc.face),1]);
%     bc.value=repmat({@(x)sin(2*pi.*x(1)).*sin(2*pi.*x(2)).*sin(2*pi.*x(3))+1},[numel(bc.face),1]);

%     perm=[kxx kxy kxz;kxy kyy kyz;kxz kyz kzz];
%     x=G.faces.centroids(:,1);y=G.faces.centroids(:,2);z=G.faces.centroids(:,3);
%     grad_u=2*pi*[cos(2*pi*x).*sin(2*pi*y).*sin(2*pi*z) sin(2*pi*x).*cos(2*pi*y).*sin(2*pi*z) sin(2*pi*x).*sin(2*pi*y).*cos(2*pi*z)]';

%     velocity=-perm*grad_u;
%     fa=dot(G.faces.normals,velocity',2);
%     hf_flux=fa(G.cells.faces(:,1));
%     hf2cn=gridCellNo(G);
%     sgn=2*(G.faces.neighbors(G.cells.faces(:,1),1)==hf2cn)-1;
%     hf_flux=sgn.*hf_flux;
%     q=accumarray(hf2cn,hf_flux);
%     src=addSource([],1:G.cells.num,q);
%     clear perm x y z grad_u velocity hf_flux sgn q hf2cn

%     disp(['Mesh refinement level ' num2str(i_con) ':'])
%     s1=OnePhaseIncompMPFA(G,TransMPFAO(G,rock,bc),src);

%     [interpFace]=findHAP(G,rock,bc);
%     hap(i_con)=interpFace.fraction;

%     interpFace=correctHAP(G,interpFace);clear coin
%     OSflux=findOSflux(G,rock,bc,interpFace);
%     u0=ones(G.cells.num,1);
%     s2=PicardNTPFA(G,bc,src,OSflux,u0,1e-7,500);
%     s3=PicardNMPFA(G,bc,src,OSflux,u0,1e-7,500);
%     picard(i_con,:)=[s2.iter s3.iter];
%     clear u0 OSflux interpFace src bc rock

%     x=G.cells.centroids(:,1);y=G.cells.centroids(:,2);z=G.cells.centroids(:,3);
%     ua=1+sin(2*pi*x).*sin(2*pi*y).*sin(2*pi*z);clear x y z
%     [ep(i_con,1),ef(i_con,1)]=L2_Norm(G,ua,fa,s1.pressure,s1.flux);
%     [ep(i_con,2),ef(i_con,2)]=L2_Norm(G,ua,fa,s2.pressure,s2.flux);
%     [ep(i_con,3),ef(i_con,3)]=L2_Norm(G,ua,fa,s3.pressure,s3.flux);
% end
% plotConvergence(ep,ef,nc,{'MPFA-O','NTPFA','NMPFA'});
% clear fa GG i_con Lx Ly myfile n 
% disp('Convergenc order of pressure solutions:')
% disp(convOrder(nc,ep,3))
% disp('Convergenc order of flux solutions:')
% disp(convOrder(nc,ef,3))
