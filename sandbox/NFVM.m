classdef NFVM < PermeabilityGradientDiscretization
    properties
        interpFace % Harmonic averaging points
        OSflux % One side fluxes
        bc % boundary conditions (own struct)
    end
    
    methods
        
        function nfvm = NFVM(model, bc, varargin) % NB BC is not optional
            
            opt = struct('myRatio', []);% tolerable ratio of harmonic average point and diameter of cell
            opt = merge_options(opt, varargin{:});
            
            % Setup nfvm members
            nfvm.bc = convertBC2FlowNTPFA(model.G, bc);
            nfvm.interpFace = nfvm.findHAP(model.G, model.rock);
            dispif(mrstVerbose, ...
                   ['fraction of faces with centroids outside convex hull: '], ...
                   num2str(nfvm.interpFace.fraction), '\n']);
            nfvm.interpFace = nfvm.correctHAP(model.G, opt.myRatio);
            nfvm.OSflux = nfvm.findOSflux(model.G, model.rock, nfvm.interpFace);
        end
        
        function v = getPermeabilityGradient(nfvm, model, state, ~)
  
            u0 = state.pressure; 
            T = nfvm.TransNTPFA(model, value(u0)); 
            v = nfvm.computeFlux(u0, T, model);
            
            % Reduce to interior
            ii = sum(model.G.faces.neighbors ~= 0, 2) == 2;
            v = -v(ii);
        end
    end
    
    methods (Access = private)
        
%         function bc2 = expandBC(nfvm, G, bc)
%             % Need to expand mrst's bc struct to get some value for _all_
%             % boundary faces. Assume the non-set bcs are homogeneous
%             % Neumann.
%             bc2.face = zeros(G.faces.num, 1);
%             bf = boundaryFaces(G);
%             bc2.face(bf) = bf;
%             bc2.type = repmat({'flux'}, [G.faces.num, 1]);
%             bc2.value = zeros(G.faces.num, 1);
%             
%             % Fill
%             bc2.type(bc.face) = bc.type;
%             bc2.value(bc.face) = bc.value;
%             
%             % Shrink
%             ii = bc2.face ~= 0;
%             bc2.face = bc2.face(ii);
%             bc2.type = bc2.type(ii);
%             bc2.value = bc2.value(ii);
%             
%             % Plot nonzero bc
%             f = zeros(G.faces.num,1);
%             f(bc2.face) = bc2.value;
%             plotGrid(G,'facealpha',0)
%             plotFaces(G,f>0);
%             xlabel 'x'
%             ylabel 'y'
%             view(3)
%         end
        
        function T = TransNTPFA(nfvm, model, u)
            dispif(mrstVerbose, 'TransNTPFA\n');
            
            G = model.G;
            % FIXME: Choose correct viscosity
            mu = model.fluid.muW(u(1));
            T=zeros(G.faces.num,2);
            for i_face=1:G.faces.num
                if(all(G.faces.neighbors(i_face,:)~=0)) % internal face
                    t1=nfvm.OSflux{i_face,1};
                    t2=nfvm.OSflux{i_face,2};
                    r1=t1(3:end-1,2)'*u(t1(3:end-1,1))+t1(end,2);
                    r2=t2(3:end-1,2)'*u(t2(3:end-1,1))+t2(end,2);
                    eps=1e-12*max(abs([t1(:,end);t2(:,end)]));
                    if(abs(r1)<=eps),r1=0;end
                    if(abs(r2)<=eps),r2=0;end
                    
                    if(abs(r1+r2)>eps)
                        mu1=r2/(r1+r2);mu2=1-mu1;
                    else
                        mu1=0.5;mu2=0.5;
                    end
                    T(i_face,1)=(mu1*t1(1,2)+mu2*t2(2,2))/mu;
                    T(i_face,2)=(mu1*t1(2,2)+mu2*t2(1,2))/mu;
                else
                    ind=find(nfvm.bc.face==i_face,1);
                    if(strcmpi(nfvm.bc.type{ind},'pressure'))
                        t1=nfvm.OSflux{i_face,1};t2=nfvm.OSflux{i_face,2};
                        t11=t1(1,2);t12=t1(2,2);
                        t22=t2(1,2);t21=t2(2,2);
                        r1=t1(3:end-1,2)'*u(t1(3:end-1,1))+t1(end,2);
                        r2=t2(end,2);
                        eps=1e-12*max(abs([t1(:,end);t2(:,end)]));
                        if(abs(r1)<=eps),r1=0;end
                        if(abs(r2)<=eps),r2=0;end
                        if(abs(r1+r2)>eps)
                            mu1=r2/(r1+r2);mu2=1-mu1;
                        else
                            mu1=0.5;mu2=0.5;
                        end
                        T(i_face,1)=mu1*t11+mu2*t21;
                        T(i_face,2)=(mu1*t12+mu2*t22)*nfvm.bc.value(ind);
                    else
                        T(i_face,2)=-G.faces.areas(i_face)*...
                            nfvm.bc.value(ind);
                    end
                end
            end
        end
        
        function flux=computeFlux(nfvm,u,T,model)
            dispif(mrstVerbose, 'computeFlux\n');
            
            G = model.G;
            
            flux=zeros(G.faces.num,1);
            flux=model.AutoDiffBackend.convertToAD(flux, u);
            ind=all(G.faces.neighbors~=0,2);
            c1=G.faces.neighbors(ind,1);c2=G.faces.neighbors(ind,2);
            flux(ind)=T(ind,1).*u(c1)-T(ind,2).*u(c2);
            c1=max(G.faces.neighbors(~ind,:),[],2);
            flux(~ind)=T(~ind,1).*u(c1)-T(~ind,2);
            ind=G.faces.neighbors(:,1)==0;
            flux(ind)=-flux(ind);
        end
        
        function interpFace=findHAP(nfvm,G,rock)
            dispif(mrstVerbose, 'findHAP\n');
            
            %find harmonic averaging points for 2D and 3D grids. Considering both
            %Dirichelt and Neumann boundary conditions
            
            %   interpFace.coords: coordinates of interpolating points
            %   interpFace.weights: interpolating weights
            %   interpFace.fraction: the fraction of cells whose centroid is
            %   outside the convex hull
            
            K=permTensor(rock,G.griddim);
            K=reshape(K',G.griddim,G.griddim,[]);
            interpFace.coords=zeros(G.faces.num,G.griddim);
            interpFace.weights=zeros(G.faces.num,2);
            interpFace.fraction=0;
            % find harmoinc averaging point--------------------------------------------
            for i_face=1:G.faces.num
                c1=G.faces.neighbors(i_face,1);
                c2=G.faces.neighbors(i_face,2);
                xf=G.faces.centroids(i_face,:)';
                if(all([c1 c2]~=0)) % internal face
                    K1=K(:,:,c1);K2=K(:,:,c2);
                    fn=G.faces.normals(i_face,:)';
                    w1=K1*fn;w2=K2*fn;
                    x1=G.cells.centroids(c1,:)';
                    x2=G.cells.centroids(c2,:)';
                    xA=x1+dot(xf-x1,fn)/dot(w1,fn)*w1;
                    xB=x2+dot(xf-x2,fn)/dot(w2,fn)*w2;
                    w1=norm(w1)/norm(xA-x1);w2=norm(w2)/norm(xB-x2);
                    interpFace.coords(i_face,:)=(w1*xA+w2*xB)'/(w1+w2);
                    interpFace.weights(i_face,1)=w1/(w1+w2);
                    interpFace.weights(i_face,2)=w2/(w1+w2);
                else
                    ind=find(nfvm.bc.face==i_face,1);
                    if(strcmpi(nfvm.bc.type{ind},'pressure'))
                        interpFace.coords(i_face,:)=xf';
                        interpFace.weights(i_face,(c2==0)+1)=nfvm.bc.value(ind);
                    else
                        c=max(c1,c2);
                        K1=K(:,:,c);
                        fn=G.faces.normals(i_face,:)';
                        w1=K1*fn;
                        x1=G.cells.centroids(c,:)';
                        xA=x1+dot(xf-x1,fn)/dot(w1,fn)*w1;
                        interpFace.coords(i_face,:)=xA';
                        a=norm(w1)/norm(x1-xA);
                        gN=nfvm.bc.value(ind);
                        interpFace.weights(i_face,(c1==0)+1)=1;
                        interpFace.weights(i_face,(c2==0)+1)=-gN/a;
                    end
                end
            end
            
            % count the number of cells whose centroid is outside the convex hull-----
            counter=zeros(G.cells.num,1);
            for i=1:G.cells.num
                xc=G.cells.centroids(i,:);
                theFaces=G.cells.faces(G.cells.facePos(i):G.cells.facePos(i+1)-1);
                hap=interpFace.coords(theFaces,:);
                ind=convhull(hap);
                switch G.griddim
                    case 2
                        xv=hap(ind,1);yv=hap(ind,2);
                        counter(i)=inpolygon(xc(1),xc(2),xv,yv);
                    case 3
                        counter(i)=nfvm.inhull(xc,hap,ind,-1e-5);
                end
            end
            interpFace.fraction=1-sum(counter)/G.cells.num;
        end
        
        function [interpFace]=correctHAP(nfvm,G,myRatio)
            dispif(mrstVerbose, 'correctHAP\n');
            
            %Correct ill-placed harmonic averaging points. If the number of input
            %arguments is 2, then the correction algorithm is applied only when some
            %cell centroids lie outside their associated convex hull; if the number of
            %input arguments is 3, then the last input argument myRatio is applied to
            %all the harmonic averaging points.
            
            %   G - Grid structure of MRST
            %   interpFace - harmonic averaging point interplation without correction
            %   myRatio - user specified ratio
            
            interpFace = nfvm.interpFace;
            
            HAP=interpFace.coords; % store the locations of the original harmonic averaging points;
            if(nargin==2 || isempty(myRatio))
                if(interpFace.fraction>0)
                    if(G.griddim==2)
                        R=0.5*G.faces.areas;
                    else
                        R=sqrt(G.faces.areas./pi);
                    end
                    flag=nfvm.isConvex(G,1:G.cells.num,interpFace);
                    while(flag)
                        mycell=flag;
                        theFaces=G.cells.faces(G.cells.facePos(mycell):G.cells.facePos(mycell+1)-1);
                        neighbors=G.faces.neighbors(theFaces,:);
                        neighbors=unique(neighbors(:));
                        neighbors(neighbors==0)=[];
                        while(flag)
                            d=interpFace.coords(theFaces,:)-G.faces.centroids(theFaces,:);
                            %d=sqrt(dot(d,d,2));
                            d = vecnorm(d,2,2);
                            [maxRatio,ind]=max(d./R(theFaces));
                            y_sigma=HAP(theFaces(ind),:)';
                            interpFace=nfvm.correctHAP_local(G,theFaces(ind),interpFace,y_sigma,0.9*maxRatio);
                            flag=nfvm.isConvex(G,mycell,interpFace);
                        end
                        flag=nfvm.isConvex(G,neighbors(1):G.cells.num,interpFace);
                    end
                end
            elseif(nargin==3)
                if(G.griddim==2)
                    R=0.5*G.faces.areas;
                else
                    R=sqrt(G.faces.areas./pi);
                end
                R=R*myRatio;
                xf=G.faces.centroids;
                hap=interpFace.coords;
                d=hap-xf;
                %d=sqrt(dot(d,d,2));
                d = vecnorm(d,2,2);
                ind=find(d>R);
                for i=1:numel(ind)
                    interpFace=correctHAP_local(G,ind(i),interpFace,HAP(ind(i),:)',myRatio);
                end
            else
                error('Wrong number of inputs')
            end
        end
        
        function flag=isConvex(nfvm,G,mycells,interpFace)
            switch G.griddim
                case 2
                    flag=0;
                    for i_cell=1:numel(mycells)
                        thecell=mycells(i_cell);
                        xc=G.cells.centroids(thecell,1);
                        yc=G.cells.centroids(thecell,2);
                        theFaces=G.cells.faces(G.cells.facePos(thecell):G.cells.facePos(thecell+1)-1);
                        hap=interpFace.coords(theFaces,:);
                        ind=convhull(hap);
                        xv=hap(ind,1);yv=hap(ind,2);
                        in=inpolygon(xc,yc,xv,yv);
                        if(~in)
                            flag=thecell;break;
                        end
                    end
                case 3
                    flag=0;
                    for i_cell=1:numel(mycells)
                        thecell=mycells(i_cell);
                        xc=G.cells.centroids(thecell,:);
                        theFaces=G.cells.faces(G.cells.facePos(thecell):G.cells.facePos(thecell+1)-1);
                        hap=interpFace.coords(theFaces,:);
                        ind=convhull(hap);
                        in=nfvm.inhull(xc,hap,ind,-1e-5);
                        %             in=inpolyhedron(ind,hap,xc);
                        if(~in),flag=thecell;break;end
                    end
            end
        end
        
        function interpFace=correctHAP_local(nfvm,G,i_face,interpFace,y_sigma,myRatio)
            % Correct harmonic averaging point for i_face based on given myRatio
            if(myRatio>0)
                if(G.griddim==2)
                    R=0.5*G.faces.areas(i_face)*myRatio;
                elseif(G.griddim==3)
                    R=myRatio*sqrt(G.faces.areas(i_face)/pi);
                end
                xm=G.faces.centroids(i_face,:)';
                interpFace.coords(i_face,:)=(xm+R*(y_sigma-xm)/norm(y_sigma-xm))';
            else
                interpFace.coords(i_face,:)=G.faces.centroids(i_face,:);
            end
        end
        
        function OSflux=findOSflux(nfvm,G,rock,interpFace)
            dispif(mrstVerbose, 'findOSflux\n');
            
            %Construct one-side fluxes for 2D and 3D grids. Considering general
            %boundary conditions, appending a constant at the last row of
            %transmissibility matrix
            % Dirichlet boundary faces are treated as zero volume cells to derive
            % nonlinear two-point flux approximation for Dirichlet boundary faces
            
            K=permTensor(rock,G.griddim);
            K=reshape(K',G.griddim,G.griddim,[]);
            OSflux=cell(G.faces.num,2);
            
            switch G.griddim
                case 2
                    for i_face=1:G.faces.num
                        if(all(G.faces.neighbors(i_face,:)~=0)) %------------------------------internal face
                            c1=G.faces.neighbors(i_face,1);
                            c2=G.faces.neighbors(i_face,2);
                            K1=K(:,:,c1);K2=K(:,:,c2);
                            w1=K1*G.faces.normals(i_face,:)';
                            w2=-K2*G.faces.normals(i_face,:)';
                            
                            [a,faceA,faceB]=nfvm.findAB(G,interpFace,c1,w1);
                            interpA=[G.faces.neighbors(faceA,:)' interpFace.weights(faceA,:)'];
                            interpB=[G.faces.neighbors(faceB,:)' interpFace.weights(faceB,:)'];
                            interpA(:,2)=-a(1)*interpA(:,2);
                            interpB(:,2)=-a(2)*interpB(:,2);
                            container=[c1;c2;interpA(:,1);interpB(:,1);0];
                            container(:,2)=[sum(a);0;interpA(:,2);interpB(:,2);0];
                            trans=nfvm.uniqueTrans(container);
                            OSflux(i_face,1)={trans};clear trans;
                            
                            [a,faceA,faceB]=nfvm.findAB(G,interpFace,c2,w2);
                            interpA=[G.faces.neighbors(faceA,:)' interpFace.weights(faceA,:)'];
                            interpB=[G.faces.neighbors(faceB,:)' interpFace.weights(faceB,:)'];
                            interpA(:,2)=-a(1)*interpA(:,2);
                            interpB(:,2)=-a(2)*interpB(:,2);
                            container=[c2;c1;interpA(:,1);interpB(:,1);0];
                            container(:,2)=[sum(a);0;interpA(:,2);interpB(:,2);0];
                            trans=nfvm.uniqueTrans(container);
                            OSflux(i_face,2)={trans};clear trans;
                        else  %--------------------------------------------boundary face
                            ind=find(nfvm.bc.face==i_face,1);
                            if(strcmpi(nfvm.bc.type{ind},'pressure'))
                                c1=max(G.faces.neighbors(i_face,:));
                                cf=G.cells.num+i_face;
                                K1=K(:,:,c1);
                                fn=G.faces.normals(i_face,:)';
                                if(c1~=G.faces.neighbors(i_face,1)),fn=-fn;end
                                w1=K1*fn;
                                [a,faceA,faceB]=nfvm.findAB(G,interpFace,c1,w1);
                                if(i_face==faceA)
                                    interpB=[G.faces.neighbors(faceB,:)' interpFace.weights(faceB,:)'];
                                    interpB(:,2)=-interpB(:,2)*a(2);
                                    container=[c1;cf;interpB(:,1);0];
                                    container(:,2)=[sum(a);-a(1);interpB(:,2);0];
                                elseif(i_face==faceB)
                                    interpA=[G.faces.neighbors(faceA,:)' interpFace.weights(faceA,:)'];
                                    interpA(:,2)=-interpA(:,2)*a(1);
                                    container=[c1;cf;interpA(:,1);0];
                                    container(:,2)=[sum(a);-a(2);interpA(:,2);0];
                                else
                                    interpA=[G.faces.neighbors(faceA,:)' interpFace.weights(faceA,:)'];
                                    interpB=[G.faces.neighbors(faceB,:)' interpFace.weights(faceB,:)'];
                                    interpA(:,2)=-a(1)*interpA(:,2);
                                    interpB(:,2)=-a(2)*interpB(:,2);
                                    container=[c1;cf;interpA(:,1);interpB(:,1);0];
                                    container(:,2)=[sum(a);0;interpA(:,2);interpB(:,2);0];
                                end
                                trans=nfvm.uniqueTrans(container);
                                OSflux(i_face,1)={trans};clear trans
                                
                                [a,xD]=nfvm.findDnode(G,c1,i_face,-w1);
                                uD=nfvm.bc.value(ind);
                                temp=[cf sum(a);c1 a(1);0 a(2)*uD];
                                OSflux(i_face,2)={temp};clear temp;
                            end
                        end
                    end
                case 3
                    for i_face=1:G.faces.num
                        if(all(G.faces.neighbors(i_face,:)~=0)) %--------------internal face
                            c1=G.faces.neighbors(i_face,1);
                            c2=G.faces.neighbors(i_face,2);

                            K1=K(:,:,c1);
                            w1=K1*G.faces.normals(i_face,:)';
                            [a,facesABC]=nfvm.findABC(G,interpFace,c1,w1);
                            OSflux(i_face,1)={nfvm.localOSflux(c1,c2,G.faces.neighbors,a,facesABC,interpFace.weights)};
                                         
                            K2=K(:,:,c2);
                            w2=-K2*G.faces.normals(i_face,:)';
                            [a,facesABC]=nfvm.findABC(G,interpFace,c2,w2);
                            OSflux(i_face,2)={nfvm.localOSflux(c2,c1,G.faces.neighbors,a,facesABC,interpFace.weights)};
                        else  %----------------------------------------------------boundary face
                            ind=find(nfvm.bc.face==i_face,1);
                            if(strcmpi(nfvm.bc.type{ind},'pressure'))
                                c1=max(G.faces.neighbors(i_face,:));
                                cf=G.cells.num+i_face;
                                K1=K(:,:,c1);fn=G.faces.normals(i_face,:)';
                                if(c1~=G.faces.neighbors(i_face,1)),fn=-fn;end
                                w1=K1*fn;
                                
                                [a,faceABC]=nfvm.findABC(G,interpFace,c1,w1);%keyboard
                                faceA=faceABC(1);
                                faceB=faceABC(2);
                                faceC=faceABC(3);
                                if(faceA==i_face)
                                    interpB=G.faces.neighbors(faceB,:)';weightB=-a(2).*interpFace.weights(faceB,:)';
                                    interpC=G.faces.neighbors(faceC,:)';weightC=-a(3).*interpFace.weights(faceC,:)';
                                    container=[c1;cf;interpB;interpC;0];
                                    container(:,2)=[sum(a);-a(1);weightB;weightC;0];
                                elseif(faceB==i_face)
                                    interpA=G.faces.neighbors(faceA,:)';weightA=-a(1).*interpFace.weights(faceA,:)';
                                    interpC=G.faces.neighbors(faceC,:)';weightC=-a(3).*interpFace.weights(faceC,:)';
                                    container=[c1;cf;interpA;interpC;0];
                                    container(:,2)=[sum(a);-a(2);weightA;weightC;0];
                                elseif(faceC==i_face)
                                    interpA=G.faces.neighbors(faceA,:)';weightA=-a(1).*interpFace.weights(faceA,:)';
                                    interpB=G.faces.neighbors(faceB,:)';weightB=-a(2).*interpFace.weights(faceB,:)';
                                    container=[c1;cf;interpA;interpB;0];
                                    container(:,2)=[sum(a);-a(3);weightA;weightB;0];
                                else
                                    interpA=G.faces.neighbors(faceA,:)';weightA=-a(1).*interpFace.weights(faceA,:)';
                                    interpB=G.faces.neighbors(faceB,:)';weightB=-a(2).*interpFace.weights(faceB,:)';
                                    interpC=G.faces.neighbors(faceC,:)';weightC=-a(3).*interpFace.weights(faceC,:)';
                                    container=[c1;cf;interpA;interpB;interpC;0];
                                    container(:,2)=[sum(a);0;weightA;weightB;weightC;0];
                                end
                                trans=nfvm.uniqueTrans(container);
                                OSflux(i_face,1)={trans};clear trans;
                                
                                [a,xA,xB]=nfvm.findDnodes(G,c1,i_face,-w1);
                                %uA=nfvm.bc.value{ind}(xA);
                                %uB=nfvm.bc.value{ind}(xB);
                                uA = nfvm.bc.value(ind);
                                uB = uA; % FIXME
                                temp=[cf sum(a);c1 a(1);0 a(2)*uA+a(3)*uB];
                                OSflux(i_face,2)={temp};clear temp;
                            end
                        end
                    end
            end
        end
        
        function trans = localOSflux(nfvm,c1,c2,neighbors,a,faces,weights)
            interpA=[neighbors(faces(1),:)' -a(1).*weights(faces(1),:)'];
            interpB=[neighbors(faces(2),:)' -a(2).*weights(faces(2),:)'];
            interpC=[neighbors(faces(3),:)' -a(3).*weights(faces(3),:)'];
            container=[c1;c2;interpA(:,1);interpB(:,1);interpC(:,1);0];
            container(:,2)=[sum(a);0;interpA(:,2);interpB(:,2);interpC(:,2);0];
            trans=nfvm.uniqueTrans(container);
        end
        
        function [a,faceA,faceB]=findAB(nfvm,G,interpFace,c,Kn)
            x1=G.cells.centroids(c,:)';
            theFaces=G.cells.faces(G.cells.facePos(c):G.cells.facePos(c+1)-1,1);
            myBases=interpFace.coords(theFaces,:);
            myBases=bsxfun(@minus,myBases,x1');
            %myNorm=sqrt(dot(myBases,myBases,2));
            myNorm=vecnorm(myBases,2,2);
            myBases=bsxfun(@rdivide,myBases,myNorm);
            Kn_norm=norm(Kn);
            Kn_unit=Kn/Kn_norm;
            myangles=bsxfun(@times,myBases,Kn_unit');
            myangles=sum(myangles,2);
            myangles=acos(myangles);
            [~,I]=sort(myangles);
            theFaces=theFaces(I);
            myBases=myBases(I,:);
            myNorm=myNorm(I);
            nf=numel(theFaces);
            flag=0;
            
            myIndex=zeros(nf*(nf-1)/2,2);
            myCoeff=myIndex;counter=1;
            for i=1:nf-1
                tA=myBases(i,:)';
                tA_norm=myNorm(i);
                for j=i+1:nf
                    tB=myBases(j,:)';
                    tB_norm=myNorm(j);
                    if(abs(det([tA tB]))>1e-9)
                        temp_a=[tA tB]\(Kn_unit);
                        temp_a(abs(temp_a)<1e-9)=0;
                        if(all(temp_a>=0))
                            if(all(temp_a<=1))
                                faceA=theFaces(i);
                                faceB=theFaces(j);
                                a=temp_a;
                                a(1)=a(1)*Kn_norm/tA_norm;
                                a(2)=a(2)*Kn_norm/tB_norm;
                                flag=1;break;
                            else
                                myIndex(counter,:)=[i,j];
                                myCoeff(counter,:)=temp_a;
                                counter=counter+1;
                            end
                        end
                    end
                end
                if(flag),break;end
            end
            if(~flag&&counter>1)
                myIndex(counter:end,:)=[];myCoeff(counter:end,:)=[];
                maxCoeff=max(myCoeff,[],2);
                [~,ind]=min(maxCoeff);
                i=myIndex(ind,1);j=myIndex(ind,2);
                a=myCoeff(ind,:);
                faceA=theFaces(i);faceB=theFaces(j);
                tA_norm=myNorm(i);
                tB_norm=myNorm(j);
                a(1)=a(1)*Kn_norm/tA_norm;
                a(2)=a(2)*Kn_norm/tB_norm;
            end
            assert(logical(exist('faceA','var')),...
                ['decomposition failed for cell ',num2str(c)]);
        end
        
        function [a,facesABC]=findABC(nfvm,G,interpFace,c,Kn)
            x1=G.cells.centroids(c,:)';
            theFaces=G.cells.faces(G.cells.facePos(c):G.cells.facePos(c+1)-1,1);
            myBases=interpFace.coords(theFaces,:);
            myBases=bsxfun(@minus,myBases,x1');
            %myNorm=sqrt(dot(myBases,myBases,2));
            myNorm=vecnorm(myBases,2,2);
            myBases=bsxfun(@rdivide,myBases,myNorm);
            Kn_norm=norm(Kn);
            Kn_unit=Kn/Kn_norm;
            myangles=bsxfun(@times,myBases,Kn_unit');
            myangles=sum(myangles,2);
            myangles=acos(myangles);
            [~,I]=sort(myangles);
            theFaces=theFaces(I);
            myBases=myBases(I,:);
            myNorm=myNorm(I);
            nf=numel(theFaces);
            flag=0;
            faces_found = false;
            
            myIndex=zeros(nf*(nf-1)*(nf-2)/6,3);
            myCoeff=myIndex;counter=1;
            for i=1:nf-2
                for j=i+1:nf-1
                    for k=j+1:nf
                        myB = myBases([i, j, k], :).';
                        if (abs(det(myB)) > 1e-9)
                            temp_a = myB \ Kn_unit;
                            temp_a(abs(temp_a)<1e-9)=0;
                            if(all(temp_a>=0))
                                if(all(temp_a<=1))
                                    facesABC = theFaces([i, j, k]);
                                    faces_found = true;
                                    a=temp_a;
                                    tABC_norm = myNorm([i, j, k]);
                                    a = a.*Kn_norm./tABC_norm;
                                    flag=1;break;
                                else
                                    myIndex(counter,:)=[i,j,k];
                                    myCoeff(counter,:)=temp_a;
                                    counter=counter+1;
                                end
                            end
                        end
                    end
                    if(flag),break;end
                end
                if(flag),break;end
            end
            if(~flag&&counter>1)
                myIndex(counter:end,:)=[];myCoeff(counter:end,:)=[];
                maxCoeff=max(myCoeff,[],2);
                [~,ind]=min(maxCoeff);
                a=myCoeff(ind,:).';
                ijk=myIndex(ind,:);
                facesABC = theFaces(ijk);
                faces_found = true;
                tABC_norm = myNorm(ijk);
                a = a.*Kn_norm./tABC_norm;
            end
            %assert(faces_found, ['decomposition failed for cell
            %',num2str(c)]);
            
            if ~faces_found
               
                figure,hold on
                plotGrid(G,c,'facealpha',0.3)
                hap=interpFace.coords(theFaces,:);
                plot3(hap(:,1),hap(:,2),hap(:,3),'.','markersize',14)
                xc=G.cells.centroids(c,:);
                plot3(xc(1),xc(2),xc(3),'o','markersize',14)
                ind=convhull(hap);
                trisurf(ind,hap(:,1),hap(:,2),hap(:,3), 'Facecolor','cyan')
                
                error(['decomposition failed for cell ', num2str(c)]);
        
            end
        end
        
        function [a,xD]=findDnode(nfvm,G,mycell,myface,Kn)
            n1=G.faces.nodes(G.faces.nodePos(myface));
            n2=G.faces.nodes(G.faces.nodePos(myface)+1);
            xn1=G.nodes.coords(n1,:)';
            xn2=G.nodes.coords(n2,:)';
            xf=G.faces.centroids(myface,:)';
            xc=G.cells.centroids(mycell,:)';
            Kn_norm=norm(Kn);
            Kn=Kn/Kn_norm;
            t_norm=norm(xc-xf);
            t=(xc-xf)/t_norm;
            t1_norm=norm(xn1-xf);
            t1=(xn1-xf)/t1_norm;
            t2_norm=norm(xn2-xf);
            t2=(xn2-xf)/t2_norm;
            temp_a=[t t1]\Kn;
            temp_a(abs(temp_a)<1e-9)=0;
            if(all(temp_a>=0))
                a=temp_a;
                a(1)=a(1)*Kn_norm/t_norm;
                a(2)=a(2)*Kn_norm/t1_norm;
                xD=xn1;
            else
                a=[t t2]\Kn;
                a(abs(a)<1e-9)=0;
                a(1)=a(1)*Kn_norm/t_norm;
                a(2)=a(2)*Kn_norm/t2_norm;
                xD=xn2;
            end
        end
        
        function [a,xA,xB]=findDnodes(nfvm,G,mycell,myface,Kn)
            mynodes=G.faces.nodes(G.faces.nodePos(myface):G.faces.nodePos(myface+1)-1);
            mynodes=[mynodes;mynodes(1)];
            xnode=G.nodes.coords(mynodes,:);
            xc=G.cells.centroids(mycell,:)';
            xf=G.faces.centroids(myface,:)';
            tc_norm=norm(xc-xf);tc=(xc-xf)/tc_norm;
            Kn_norm=norm(Kn);Kn=Kn/Kn_norm;
            for i=1:numel(mynodes)-1
                xA=xnode(i,:)';
                xB=xnode(i+1,:)';
                tA_norm=norm(xA-xf);tA=(xA-xf)/tA_norm;
                tB_norm=norm(xB-xf);tB=(xB-xf)/tB_norm;
                a=[tc tA tB]\Kn;
                a(abs(a)<1e-9)=0;
                if(all(a>=0))
                    a(1)=a(1)*Kn_norm/tc_norm;
                    a(2)=a(2)*Kn_norm/tA_norm;
                    a(3)=a(3)*Kn_norm/tB_norm;
                    break;
                end
            end
        end
        
        function [trans]=uniqueTrans(nfvm,container)
            [trans,~,subs]=unique(container(:,1),'rows','stable');
            trans(:,2)=accumarray(subs,container(:,2));
            trans(3:end,:)=sortrows(trans(3:end,:),-1);
            trans(2:end,2)=-trans(2:end,2);
        end
        
        function in = inhull(nfvm,testpts,xyz,tess,tol)
            % inhull: tests if a set of points are inside a convex hull
            % usage: in = inhull(testpts,xyz)
            % usage: in = inhull(testpts,xyz,tess)
            % usage: in = inhull(testpts,xyz,tess,tol)
            %
            % arguments: (input)
            %  testpts - nxp array to test, n data points, in p dimensions
            %       If you have many points to test, it is most efficient to
            %       call this function once with the entire set.
            %
            %  xyz - mxp array of vertices of the convex hull, as used by
            %       convhulln.
            %
            %  tess - tessellation (or triangulation) generated by convhulln
            %       If tess is left empty or not supplied, then it will be
            %       generated.
            %
            %  tol - (OPTIONAL) tolerance on the tests for inclusion in the
            %       convex hull. You can think of tol as the distance a point
            %       may possibly lie outside the hull, and still be perceived
            %       as on the surface of the hull. Because of numerical slop
            %       nothing can ever be done exactly here. I might guess a
            %       semi-intelligent value of tol to be
            %
            %         tol = 1.e-13*mean(abs(xyz(:)))
            %
            %       In higher dimensions, the numerical issues of floating
            %       point arithmetic will probably suggest a larger value
            %       of tol.
            %
            %       DEFAULT: tol = 0
            %
            % arguments: (output)
            %  in  - nx1 logical vector
            %        in(i) == 1 --> the i'th point was inside the convex hull.
            %
            % Example usage: The first point should be inside, the second out
            %
            %  xy = randn(20,2);
            %  tess = convhulln(xy);
            %  testpoints = [ 0 0; 10 10];
            %  in = inhull(testpoints,xy,tess)
            %
            % in =
            %      1
            %      0
            %
            % A non-zero count of the number of degenerate simplexes in the hull
            % will generate a warning (in 4 or more dimensions.) This warning
            % may be disabled off with the command:
            %
            %   warning('off','inhull:degeneracy')
            %
            % See also: convhull, convhulln, delaunay, delaunayn, tsearch, tsearchn
            %
            % Author: John D'Errico
            % e-mail: woodchips@rochester.rr.com
            % Release: 3.0
            % Release date: 10/26/06
            
            % get array sizes
            % m points, p dimensions
            p = size(xyz,2);
            [n,c] = size(testpts);
            if p ~= c
                error 'testpts and xyz must have the same number of columns'
            end
            if p < 2
                error 'Points must lie in at least a 2-d space.'
            end
            
            % was the convex hull supplied?
            if (nargin<3) || isempty(tess)
                tess = convhulln(xyz);
            end
            [nt,c] = size(tess);
            if c ~= p
                error 'tess array is incompatible with a dimension p space'
            end
            
            % was tol supplied?
            if (nargin<4) || isempty(tol)
                tol = 0;
            end
            
            % build normal vectors
            switch p
                case 2
                    % really simple for 2-d
                    nrmls = (xyz(tess(:,1),:) - xyz(tess(:,2),:)) * [0 1;-1 0];
                    
                    % Any degenerate edges?
                    del = sqrt(sum(nrmls.^2,2));
                    degenflag = (del<(max(del)*10*eps));
                    if sum(degenflag)>0
                        warning('inhull:degeneracy',[num2str(sum(degenflag)), ...
                            ' degenerate edges identified in the convex hull'])
                        
                        % we need to delete those degenerate normal vectors
                        nrmls(degenflag,:) = [];
                        nt = size(nrmls,1);
                    end
                case 3
                    % use vectorized cross product for 3-d
                    ab = xyz(tess(:,1),:) - xyz(tess(:,2),:);
                    ac = xyz(tess(:,1),:) - xyz(tess(:,3),:);
                    nrmls = cross(ab,ac,2);
                    degenflag = false(nt,1);
                otherwise
                    % slightly more work in higher dimensions,
                    nrmls = zeros(nt,p);
                    degenflag = false(nt,1);
                    for i = 1:nt
                        % just in case of a degeneracy
                        % Note that bsxfun COULD be used in this line, but I have chosen to
                        % not do so to maintain compatibility. This code is still used by
                        % users of older releases.
                        %  nullsp = null(bsxfun(@minus,xyz(tess(i,2:end),:),xyz(tess(i,1),:)))';
                        nullsp = null(xyz(tess(i,2:end),:) - repmat(xyz(tess(i,1),:),p-1,1))';
                        if size(nullsp,1)>1
                            degenflag(i) = true;
                            nrmls(i,:) = NaN;
                        else
                            nrmls(i,:) = nullsp;
                        end
                    end
                    if sum(degenflag)>0
                        warning('inhull:degeneracy',[num2str(sum(degenflag)), ...
                            ' degenerate simplexes identified in the convex hull'])
                        
                        % we need to delete those degenerate normal vectors
                        nrmls(degenflag,:) = [];
                        nt = size(nrmls,1);
                    end
            end
            
            % scale normal vectors to unit length
            nrmllen = sqrt(sum(nrmls.^2,2));
            % again, bsxfun COULD be employed here...
            %  nrmls = bsxfun(@times,nrmls,1./nrmllen);
            nrmls = nrmls.*repmat(1./nrmllen,1,p);
            
            % center point in the hull
            center = mean(xyz,1);
            
            % any point in the plane of each simplex in the convex hull
            a = xyz(tess(~degenflag,1),:);
            
            % ensure the normals are pointing inwards
            % this line too could employ bsxfun...
            %  dp = sum(bsxfun(@minus,center,a).*nrmls,2);
            dp = sum((repmat(center,nt,1) - a).*nrmls,2);
            k = dp<0;
            nrmls(k,:) = -nrmls(k,:);
            
            % We want to test if:  dot((x - a),N) >= 0
            % If so for all faces of the hull, then x is inside
            % the hull. Change this to dot(x,N) >= dot(a,N)
            aN = sum(nrmls.*a,2);
            
            % test, be careful in case there are many points
            in = false(n,1);
            
            % if n is too large, we need to worry about the
            % dot product grabbing huge chunks of memory.
            memblock = 1e6;
            blocks = max(1,floor(n/(memblock/nt)));
            aNr = repmat(aN,1,length(1:blocks:n));
            for i = 1:blocks
                j = i:blocks:n;
                if size(aNr,2) ~= length(j),
                    aNr = repmat(aN,1,length(j));
                end
                in(j) = all((nrmls*testpts(j,:)' - aNr) >= -tol,1)';
            end
        end
        
    end
end
