function [ T ] = TransMPFAO(Grid,rock,bc,varargin)
%Compute face transmissiblity using MPFA-O. Applicable to general 2-D
%polygonal grids and 3-D polyhedral grids (when MPFA-O is applicable)
%   Grid - Grid structure of MRST
%   rock - Rock structure of MRST
%   bc - Boundary condition structure.
%   T - Transmissibility cell array, one 2-D matrix for each face

K=permTensor(rock,Grid.griddim);
K=reshape(K',Grid.griddim,Grid.griddim,[]);
nn=Grid.nodes.num;
nf=Grid.faces.num;
if(nargin>3)
    param=varargin{:};
else
    param=0;
end
switch Grid.griddim
    case 2   % ------------------------------------------2D polygonal grids
        TT=cell(nn,1);
        for i_node=1:nn
            xO=Grid.nodes.coords(i_node,:)';
            [Trans,counter]=nodeInfo(Grid,i_node);
            
            if(all(Trans(:,2)~=0)) %-------------------internal node
                x_cell=Grid.cells.centroids(Trans(:,2),:);
                xm=Grid.faces.centroids(Trans(:,1),:);
                x_face=xm+param.*bsxfun(@minus,xO',xm);
                face_normal=0.5*Grid.faces.normals(Trans(:,1),:);
                perm=K(:,:,Trans(:,2));
                
                M=zeros(2,3,counter);
                for i=1:counter
                    M(:,:,i)=GradTri(x_cell(i,:)',x_face(i,:)',...
                        x_face(stepBack(i,counter),:)');
                end
                
                G=zeros(2,3,counter);
                for i=1:counter
                    G(:,:,i)=-[face_normal(i,:)*perm(:,:,i)*M(:,:,i);...
                        face_normal(stepBack(i,counter),:)...
                        *perm(:,:,i)*M(:,:,i)];
                end
                A=zeros(counter,counter);
                B=A;
                C=A;
                D=A;
                for i=1:counter
                    A(i,i)=G(1,2,i)-G(2,3,stepFor(i,counter));
                    A(i,stepFor(i,counter))=-G(2,2,stepFor(i,counter));
                    A(i,stepBack(i,counter))=G(1,3,i);
                    B(i,i)=-G(1,1,i);
                    B(i,stepFor(i,counter))=G(2,1,stepFor(i,counter));
                    C(i,i)=G(1,2,i);
                    C(i,stepBack(i,counter))=G(1,3,i);
                    D(i,i)=G(1,1,i);
                end
                Trans(:,3:end)=C*A^-1*B+D;
                TT{i_node}=Trans;
            else  %---------------------------------------Boundary node
                x_cell=Grid.cells.centroids(Trans(1:end-1,2),:);
                xm=Grid.faces.centroids(Trans(:,1),:);
                x_face=xm+param.*bsxfun(@minus,xO',xm);
                face_normal=0.5*Grid.faces.normals(Trans(:,1),:);
                perm=K(:,:,Trans(1:end-1,2));
                
                M=zeros(2,3,counter-1);
                for i=1:counter-1
                    M(:,:,i)=GradTri(x_cell(i,:)',x_face(i,:)',...
                        x_face(stepBack(i,counter),:)');
                end
                
                G=zeros(2,3,counter-1);
                for i=1:counter-1
                    G(:,:,i)=-[face_normal(i,:)*perm(:,:,i)*M(:,:,i);...
                        face_normal(stepBack(i,counter),:)*...
                        perm(:,:,i)*M(:,:,i)];
                end
                
                A=zeros(counter,counter);
                B=A;
                C=A;
                D=A;
                for i=1:counter-2
                    A(i,i)=G(1,2,i)-G(2,3,i+1);
                    A(i,i+1)=-G(2,2,i+1);
                    A(i,stepBack(i,counter))=G(1,3,i);
                    B(i,i)=-G(1,1,i);
                    B(i,i+1)=G(2,1,i+1);
                    C(i,i)=G(1,2,i);
                    C(i,stepBack(i,counter))=G(1,3,i);
                    D(i,i)=G(1,1,i);
                end
                
                indexC=find(bc.face==Trans(end-1,1),1);
                indexD=find(bc.face==Trans(end,1),1);
                dataC=bc.value{indexC};
                dataD=bc.value{indexD};
                
                if(strcmpi(bc.type{indexC},'pressure'))
                    A(end-1,end-1)=1;
                    B(end-1,end)=dataC(x_face(end-1,:));
                else
                    A(counter-1,stepBack(counter-1,counter))=G(1,3,end);
                    A(counter-1,counter-1)=G(1,2,end);
                    B(counter-1,counter-1)=-G(1,1,end);
                    gN=0.5*dataC(xm(end-1,:))*Grid.faces.areas(Trans(end-1,1));
                    if(Grid.faces.neighbors(Trans(end-1,1),1)==0)
                        gN=-gN;
                    end
                    B(counter-1,end)=gN;
                end
                
                if(strcmp(bc.type{indexD},'pressure'))
                    A(end,end)=1;
                    B(end,end)=dataD(x_face(end,:));
                else
                    A(end,1)=G(2,2,1);
                    A(end,end)=G(2,3,1);
                    B(end,1)=-G(2,1,1);
                    gN=0.5*dataD(xm(end,:))*Grid.faces.areas(Trans(end,1));
                    if(Grid.faces.neighbors(Trans(end,1),1)==0)
                        gN=-gN;
                    end
                    B(end,end)=gN;
                end
                
                C(counter-1,counter-1)=G(1,2,end);
                C(counter-1,stepBack(counter-1,counter))=G(1,3,end);
                C(counter,counter)=G(2,3,1);
                C(counter,1)=G(2,2,1);
                D(end-1,end-1)=G(1,1,end);
                D(end,1)=G(2,1,1);
                Trans(:,3:end)=C*A^-1*B+D;
                TT{i_node}=Trans;
            end
        end
        
        T=cell(nf,1);
        for i_face=1:nf
            n1=Grid.faces.nodes(Grid.faces.nodePos(i_face));
            n2=Grid.faces.nodes(Grid.faces.nodePos(i_face)+1);
            T1=TT{n1};
            T2=TT{n2};
            ind1=find(i_face==T1(:,1),1);
            ind2=find(i_face==T2(:,1),1);
            container=[T1(:,2) T1(ind1,3:end)';T2(:,2) T2(ind2,3:end)'];
            [temp,~,subs]=unique(container(:,1),'rows','stable');
            temp(:,2)=accumarray(subs,container(:,2));
            temp=sortrows(temp,-1);
            T{i_face}=temp;
        end
    case 3   % -----------------------------------------3D polyhedral grids
        TT=repmat(struct('faces',[],'cells',[],'Trans',[]),nn,1);
        for i_node=1:nn
            I=find(Grid.faces.nodes==i_node);
            counter=length(I);
            myfaces=zeros(counter,1);
            for i=1:counter
                myfaces(i)=find(Grid.faces.nodePos(1:end-1)<=I(i)&...
                    Grid.faces.nodePos(2:end)>I(i));
            end
            clear i I counter
            mycells=Grid.faces.neighbors(myfaces,:);
            mycells=unique(mycells(:));
            mycells=sortrows(mycells,-1);
            TT(i_node).faces=myfaces;
            TT(i_node).cells=mycells;
            
            cellInfo=repmat(struct('faces',[],'G',[]),nnz(mycells),1);
            for i=1:nnz(mycells)
                c=mycells(i);
                cfaces=Grid.cells.faces(Grid.cells.facePos(c):...
                    Grid.cells.facePos(c+1)-1);
                cfaces=intersect(myfaces,cfaces);
                assert(numel(cfaces)==3);
                cellInfo(i).faces=cfaces;
                fA=cfaces(1);fB=cfaces(2);fC=cfaces(3);
                nfA=Grid.faces.nodePos(fA+1)-Grid.faces.nodePos(fA);
                nfB=Grid.faces.nodePos(fB+1)-Grid.faces.nodePos(fB);
                nfC=Grid.faces.nodePos(fC+1)-Grid.faces.nodePos(fC);
                xi=Grid.cells.centroids(c,:)';
                xA=Grid.faces.centroids(fA,:)';
                xB=Grid.faces.centroids(fB,:)';
                xC=Grid.faces.centroids(fC,:)';
                fnA=Grid.faces.normals(fA,:)'/nfA;
                fnB=Grid.faces.normals(fB,:)'/nfB;
                fnC=Grid.faces.normals(fC,:)'/nfC;
                [dNi,dNA,dNB,dNC]=GradTet(xi,xA,xB,xC);
                perm=K(:,:,c);
                cellInfo(i).G=[-fnA';-fnB';-fnC']*perm*[dNi dNA dNB dNC];
            end
            
            A=zeros(numel(myfaces));
            B=zeros(numel(myfaces),numel(mycells));
            C=zeros(numel(myfaces));
            D=zeros(numel(myfaces),numel(mycells));
            for i=1:numel(myfaces)
                fi=myfaces(i);
                if(all(Grid.faces.neighbors(fi,:)~=0))
                    c1=Grid.faces.neighbors(fi,1);
                    c2=Grid.faces.neighbors(fi,2);
                    ind_c1=find(c1==mycells,1);
                    ind_c2=find(c2==mycells,1);
                    f1A=cellInfo(ind_c1).faces(1);
                    ind_f1A=find(f1A==myfaces,1);
                    f1B=cellInfo(ind_c1).faces(2);
                    ind_f1B=find(f1B==myfaces,1);
                    f1C=cellInfo(ind_c1).faces(3);
                    ind_f1C=find(f1C==myfaces,1);
                    f2A=cellInfo(ind_c2).faces(1);
                    ind_f2A=find(f2A==myfaces,1);
                    f2B=cellInfo(ind_c2).faces(2);
                    ind_f2B=find(f2B==myfaces,1);
                    f2C=cellInfo(ind_c2).faces(3);
                    ind_f2C=find(f2C==myfaces,1);
                    G=[cellInfo(ind_c1).G(cellInfo(ind_c1).faces==fi,:);...
                        cellInfo(ind_c2).G(cellInfo(ind_c2).faces==fi,:)];
                    A(i,[ind_f1A ind_f1B ind_f1C])=...
                        A(i,[ind_f1A ind_f1B ind_f1C])+G(1,2:4);
                    A(i,[ind_f2A ind_f2B ind_f2C])=...
                        A(i,[ind_f2A ind_f2B ind_f2C])-G(2,2:4);
                    B(i,[ind_c1 ind_c2])=[-G(1,1) G(2,1)];
                    C(i,[ind_f1A ind_f1B ind_f1C])=G(1,2:4);
                    D(i,ind_c1)=G(1,1);
                else
                    c1=max(Grid.faces.neighbors(fi,:));
                    ind_c1=find(c1==mycells,1);
                    f1A=cellInfo(ind_c1).faces(1);
                    ind_f1A=find(f1A==myfaces,1);
                    f1B=cellInfo(ind_c1).faces(2);
                    ind_f1B=find(f1B==myfaces,1);
                    f1C=cellInfo(ind_c1).faces(3);
                    ind_f1C=find(f1C==myfaces,1);
                    G=cellInfo(ind_c1).G(cellInfo(ind_c1).faces==fi,:);
                    
                    ind=find(bc.face==fi,1);
                    data=bc.value{ind};
                    if(strcmpi(bc.type{ind},'pressure'))
                        xf=Grid.faces.centroids(fi,:)';
                        A(i,i)=1;
                        B(i,end)=data(xf);
                        C(i,[ind_f1A ind_f1B ind_f1C])=G(2:4);
                        D(i,ind_c1)=G(1);
                    else
                        A(i,[ind_f1A ind_f1B ind_f1C])=G(2:4);
                        B(i,ind_c1)=-G(1);
                        temp=Grid.faces.nodePos(fi+1)-Grid.faces.nodePos(fi);
                        gN=data(Grid.faces.centroids(fi,:))*Grid.faces.areas(fi);
                        gN=gN/temp;
                        if(Grid.faces.neighbors(fi,1)==0)
                            gN=-gN;
                        end
                        B(i,end)=gN;
                        D(i,end)=gN;
                    end
                end
            end
            
            TT(i_node).Trans=C*A^-1*B+D;
        end
        
        T=cell(nf,1);
        for i_face=1:nf
            mynodes=Grid.faces.nodes(Grid.faces.nodePos(i_face):...
                Grid.faces.nodePos(i_face+1)-1);
            container=cell(numel(mynodes),1);
            for i_node=1:numel(mynodes)
                t=TT(mynodes(i_node)).Trans;
                container(i_node)={[TT(mynodes(i_node)).cells ...
                    t(i_face==TT(mynodes(i_node)).faces,:)']};
            end
            container=cell2mat(container);
            [temp,~,subs]=unique(container(:,1),'rows','stable');
            temp(:,2)=accumarray(subs,container(:,2));
            temp=sortrows(temp,-1);
            T{i_face}=temp;
        end
end
end

%% Auxilliary functions
function [M,A]=GradTri(xi,xj,xk)
bi=xj(2)-xk(2);bj=xk(2)-xi(2);bk=xi(2)-xj(2);
ci=xk(1)-xj(1);cj=xi(1)-xk(1);ck=xj(1)-xi(1);
A=abs(det([1 xi';1 xj';1 xk']));
M=[bi bj bk;ci cj ck]/A;
end
function [dNi,dNj,dNk,dNl]=GradTet(xi,xj,xk,xl)
V=abs(det([xi-xj,xi-xk,xi-xl]));
dNi=((2*(dot(cross(xk-xj,xl-xj),xi-xj)>0)-1)*cross(xk-xj,xl-xj))/V;
dNj=((2*(dot(cross(xk-xi,xl-xi),xj-xi)>0)-1)*cross(xk-xi,xl-xi))/V;
dNk=((2*(dot(cross(xj-xi,xl-xi),xk-xi)>0)-1)*cross(xj-xi,xl-xi))/V;
dNl=((2*(dot(cross(xj-xi,xk-xi),xl-xi)>0)-1)*cross(xj-xi,xk-xi))/V;
end
function y=stepBack(i,n)
if(i==1)
    y=n;
else
    y=i-1;
end
end
function y=stepFor(i,n)
if(i<n)
    y=i+1;
else
    y=1;
end
end
function [Trans,counter]=nodeInfo(Grid,i_node)
% find faces and cells connected to i_node and put them in
% counter-clockwise order for 2-D polygonal grids
R=[0 -1;1 0];
temp(:,1)=ceil(find(Grid.faces.nodes==i_node)*0.5);
temp(:,2:3)=Grid.faces.neighbors(temp(:,1),:);
counter=size(temp,1);
myfaces=temp(:,1);

Trans=zeros(counter,counter+2);
xO=Grid.nodes.coords(i_node,:)';
if(all([temp(:,2);temp(:,3)]~=0)) %--------------------------internal node
    Trans(1,1)=temp(1,1);
    c1=temp(1,2);c2=temp(1,3);
    clear temp;
    xm=Grid.faces.centroids(Trans(1,1),:)';
    x1=Grid.cells.centroids(c1,:)';
    x2=Grid.cells.centroids(c2,:)';
    if((xO-xm)'*R*(x2-x1)>0)
        Trans(1,2)=c1;Trans(2,2)=c2;
    else
        Trans(1,2)=c2;Trans(2,2)=c1;
    end
    for i=2:counter
        c=Trans(i,2);
        cfaces=Grid.cells.faces(Grid.cells.facePos(c):...
            Grid.cells.facePos(c+1)-1);
        cfaces=intersect(cfaces,myfaces);
        cfaces(cfaces==Trans(i-1,1))=[];
        Trans(i,1)=cfaces;
        if(i<counter)
            thecells=Grid.faces.neighbors(Trans(i,1),:);
            thecells(thecells==c)=[];
            Trans(i+1,2)=thecells;
        end
    end
else  %------------------------------------------------------boundary node
    ind=any(temp==0,2);
    temp(~ind,:)=[];
    xm=Grid.faces.centroids(temp(1,1),:)';
    c=max(temp(1,2:3));
    xc=Grid.cells.centroids(c,:)';
    if((xc-xO)'*R*(xm-xO)>0)
        Trans(end,1)=temp(1,1);
        Trans(end-1,1)=temp(2,1);
        Trans(1,2)=max(temp(1,2:3));
        Trans(end-1,2)=max(temp(2,2:3));
    else
        Trans(end,1)=temp(2,1);
        Trans(end-1,1)=temp(1,1);
        Trans(1,2)=max(temp(2,2:3));
        Trans(end-1,2)=max(temp(1,2:3));
    end
    for i=1:counter-2
        c=Trans(i,2);
        cfaces=Grid.cells.faces(Grid.cells.facePos(c):...
            Grid.cells.facePos(c+1)-1);
        cfaces=intersect(cfaces,myfaces);
        if(any(cfaces(1)==Trans(:,1)))
            Trans(i,1)=cfaces(2);
        else
            Trans(i,1)=cfaces(1);
        end
        thecells=Grid.faces.neighbors(Trans(i,1),:);
        thecells(thecells==Trans(i,2))=[];
        Trans(i+1,2)=thecells;
    end
end
end
