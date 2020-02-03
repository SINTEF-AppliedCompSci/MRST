function OSflux=findOSflux(G,rock,bc,interpFace)
%Construct one-side fluxes for 2D and 3D grids. Considering general
%boundary conditions, appending a constant at the last row of
%transmissibility matrix
% Dirichlet boundary faces are treated as zero volume cells to derive
% nonlinear two-point flux approximation for Dirichlet boundary faces
dispif(mrstVerbose, 'findOSflux\n');

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
                
                [a,faceA,faceB]=findAB(G,interpFace,c1,w1);
                interpA=[G.faces.neighbors(faceA,:)' interpFace.weights(faceA,:)'];
                interpB=[G.faces.neighbors(faceB,:)' interpFace.weights(faceB,:)'];
                interpA(:,2)=-a(1)*interpA(:,2);
                interpB(:,2)=-a(2)*interpB(:,2);
                container=[c1;c2;interpA(:,1);interpB(:,1);0];
                container(:,2)=[sum(a);0;interpA(:,2);interpB(:,2);0];
                trans=uniqueTrans(container);
                OSflux(i_face,1)={trans};clear trans;

                [a,faceA,faceB]=findAB(G,interpFace,c2,w2);
                interpA=[G.faces.neighbors(faceA,:)' interpFace.weights(faceA,:)'];
                interpB=[G.faces.neighbors(faceB,:)' interpFace.weights(faceB,:)'];
                interpA(:,2)=-a(1)*interpA(:,2);
                interpB(:,2)=-a(2)*interpB(:,2);
                container=[c2;c1;interpA(:,1);interpB(:,1);0];
                container(:,2)=[sum(a);0;interpA(:,2);interpB(:,2);0];
                trans=uniqueTrans(container);
                OSflux(i_face,2)={trans};clear trans;
            else  %--------------------------------------------boudary face
                ind=find(bc.face==i_face,1);
                if(strcmpi(bc.type{ind},'pressure'))
                    c1=max(G.faces.neighbors(i_face,:));
                    cf=G.cells.num+i_face;
                    K1=K(:,:,c1);
                    fn=G.faces.normals(i_face,:)';
                    if(c1~=G.faces.neighbors(i_face,1)),fn=-fn;end
                    w1=K1*fn;
                    [a,faceA,faceB]=findAB(G,interpFace,c1,w1);
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
                    trans=uniqueTrans(container);
                    OSflux(i_face,1)={trans};clear trans
                    
                    [a,xD]=findDnode(G,c1,i_face,-w1);
                    uD=bc.value{ind}(xD);
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
                K1=K(:,:,c1);K2=K(:,:,c2);
                w1=K1*G.faces.normals(i_face,:)';
                w2=-K2*G.faces.normals(i_face,:)';
                
                [a,faceA,faceB,faceC]=findABC(G,interpFace,c1,w1);
                interpA=[G.faces.neighbors(faceA,:)' -a(1).*interpFace.weights(faceA,:)'];
                interpB=[G.faces.neighbors(faceB,:)' -a(2).*interpFace.weights(faceB,:)'];
                interpC=[G.faces.neighbors(faceC,:)' -a(3).*interpFace.weights(faceC,:)'];
                container=[c1;c2;interpA(:,1);interpB(:,1);interpC(:,1);0];
                container(:,2)=[sum(a);0;interpA(:,2);interpB(:,2);interpC(:,2);0];
                trans=uniqueTrans(container);
                OSflux(i_face,1)={trans};clear trans;

                [a,faceA,faceB,faceC]=findABC(G,interpFace,c2,w2);
                interpA=[G.faces.neighbors(faceA,:)' -a(1).*interpFace.weights(faceA,:)'];
                interpB=[G.faces.neighbors(faceB,:)' -a(2).*interpFace.weights(faceB,:)'];
                interpC=[G.faces.neighbors(faceC,:)' -a(3).*interpFace.weights(faceC,:)'];
                container=[c2;c1;interpA(:,1);interpB(:,1);interpC(:,1);0];
                container(:,2)=[sum(a);0;interpA(:,2);interpB(:,2);interpC(:,2);0];
                trans=uniqueTrans(container);
                OSflux(i_face,2)={trans};clear trans;
            else  %----------------------------------------------------boudary face
                ind=find(bc.face==i_face,1);
                if(strcmpi(bc.type{ind},'pressure'))
                    c1=max(G.faces.neighbors(i_face,:));
                    cf=G.cells.num+i_face;
                    K1=K(:,:,c1);fn=G.faces.normals(i_face,:)';
                    if(c1~=G.faces.neighbors(i_face,1)),fn=-fn;end
                    w1=K1*fn;
                    
                    [a,faceA,faceB,faceC]=findABC(G,interpFace,c1,w1);
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
                    trans=uniqueTrans(container);
                    OSflux(i_face,1)={trans};clear trans;

                    [a,xA,xB]=findDnodes(G,c1,i_face,-w1);
                    uA=bc.value{ind}(xA);
                    uB=bc.value{ind}(xB);
                    temp=[cf sum(a);c1 a(1);0 a(2)*uA+a(3)*uB];
                    OSflux(i_face,2)={temp};clear temp;
                end
            end
        end
end
end

function [a,faceA,faceB]=findAB(G,interpFace,c,Kn)
x1=G.cells.centroids(c,:)';
theFaces=G.cells.faces(G.cells.facePos(c):G.cells.facePos(c+1)-1,1);
myBases=interpFace.coords(theFaces,:);
myBases=bsxfun(@minus,myBases,x1');
myNorm=sqrt(dot(myBases,myBases,2));
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

function [a,faceA,faceB,faceC]=findABC(G,interpFace,c,Kn)
x1=G.cells.centroids(c,:)';
theFaces=G.cells.faces(G.cells.facePos(c):G.cells.facePos(c+1)-1,1);
myBases=interpFace.coords(theFaces,:);
myBases=bsxfun(@minus,myBases,x1');
myNorm=sqrt(dot(myBases,myBases,2));
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

% if mrstDebug
%     figure, hold on
%     plotGrid(G,c,'facealpha',0.1);
%     [Gc,~,~,gn] = extractSubgrid(G,c);
%     xmin = min(G.nodes.coords(gn,:));
%     xmax = max(G.nodes.coords(gn,:));
%     h = xmax-xmin;
%     plot3(x1(1),x1(2),x1(3),'.','markersize',15);
%     for i = 1:numel(theFaces)
%         X = [x1(1); x1(1)+myBases(i,1)*h(1)];
%         Y = [x1(2); x1(2)+myBases(i,2)*h(2)];
%         Z = [x1(3); x1(3)+myBases(i,3)*h(3)];
%         line(X,Y,Z);
%     end
%     keyboard
% end
%faces_found = false;

myIndex=zeros(nf*(nf-1)*(nf-2)/6,3);
myCoeff=myIndex;counter=1;
for i=1:nf-2
    tA=myBases(i,:)';
    tA_norm=myNorm(i);
    for j=i+1:nf-1
        tB=myBases(j,:)';
        tB_norm=myNorm(j);
        for k=j+1:nf
            tC=myBases(k,:)';
            tC_norm=myNorm(k);
            %det([tA tB tC])
            if(abs(det([tA tB tC]))>1e-9)
                temp_a=[tA tB tC]\(Kn_unit);
                temp_a(abs(temp_a)<1e-9)=0;
                %temp_a
                if(all(temp_a>=0))
                    if(all(temp_a<=1))
                        faceA=theFaces(i);
                        faceB=theFaces(j);
                        faceC=theFaces(k);
                        a=temp_a;
                        a(1)=a(1)*Kn_norm/tA_norm;
                        a(2)=a(2)*Kn_norm/tB_norm;
                        a(3)=a(3)*Kn_norm/tC_norm;
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
    i=myIndex(ind,1);j=myIndex(ind,2);k=myIndex(ind,3);
    a=myCoeff(ind,:);
    faceA=theFaces(i);faceB=theFaces(j);faceC=theFaces(k);
    tA_norm=myNorm(i);
    tB_norm=myNorm(j);
    tC_norm=myNorm(k);
    a(1)=a(1)*Kn_norm/tA_norm;
    a(2)=a(2)*Kn_norm/tB_norm;
    a(3)=a(3)*Kn_norm/tC_norm;
        
    %keyboard
    
    %faces_found = true;
end

% if ~faces_found
%     
%     figure,hold on
%     plotGrid(G,c,'facealpha',0.3)
%     hap=interpFace.coords(theFaces,:);
%     plot3(hap(:,1),hap(:,2),hap(:,3),'.','markersize',14)
%     xc=G.cells.centroids(c,:);
%     plot3(xc(1),xc(2),xc(3),'o','markersize',14)
%     ind=convhull(hap);
%     trisurf(ind,hap(:,1),hap(:,2),hap(:,3), 'Facecolor','cyan')
%     
%     figure
%     plotGrid(G, 'facealpha', 0.1);
%     hold on
%     plotGrid(G, c)
% 
%     keyboard
% 
%     disp(['decomposition failed for cell ', num2str(c)]);
%     %error(['decomposition failed for cell ', num2str(c)]);
%     
% end

if ~exist('faceA','var') %|| nf == 4
    
    figure, hold on
    plotGrid(G, 'facealpha', 0.1);
    plotGrid(G, c)

    figure,hold on
    plotGrid(G,c,'facealpha',0.3)
    hap=interpFace.coords(theFaces,:);
    plot3(hap(:,1),hap(:,2),hap(:,3),'.','markersize',14)
    xc=G.cells.centroids(c,:);
    plot3(xc(1),xc(2),xc(3),'ko','markersize',14)
    ind=convhull(hap);
    trisurf(ind,hap(:,1),hap(:,2),hap(:,3), 'Facecolor','cyan','facealpha',0.5)
    in_convex_hull=max_inhull(xc,hap,ind,1e-5)
    
    % Plot Kn. 
    Gc = extractSubgrid(G,c);
    xmin = min(Gc.nodes.coords);
    xmax = max(Gc.nodes.coords);
    h = xmax-xmin;
    vec = [xc; xc+Kn_unit'.*h];
    line(vec(:,1),vec(:,2),vec(:,3))
    
    keyboard
    assert(logical(exist('faceA','var')),...
       ['decomposition failed for cell ',num2str(c)]);
   
%    % Take three faces with largest alignment with Kn. If we end up here,
%    % they're typically not all positive.
%    disp(['fix faceA,B,C ', num2str(c), ' ', num2str(G.cells.num)])
%    dotp = sum(myBases.*Kn', 2);
%    [dotp, ii] = sort(dotp);
%    faceA = theFaces(ii(1));
%    faceB = theFaces(ii(2));
%    faceC = theFaces(ii(3));
%    tA = myBases(ii(1),:)';
%    tB = myBases(ii(2),:)';
%    tC = myBases(ii(3),:)';
%    a=[tA tB tC]\(Kn_unit);
%    a(1)=a(1)*Kn_norm/norm(tA);
%    a(2)=a(2)*Kn_norm/norm(tB);
%    a(3)=a(3)*Kn_norm/norm(tC);
%    %keyboard
    
end

end

function [a,xD]=findDnode(G,mycell,myface,Kn)
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

function [a,xA,xB]=findDnodes(G,mycell,myface,Kn)
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

function [trans]=uniqueTrans(container)
[trans,~,subs]=unique(container(:,1),'rows','stable');
trans(:,2)=accumarray(subs,container(:,2));
trans(3:end,:)=sortrows(trans(3:end,:),-1);
trans(2:end,2)=-trans(2:end,2);
end


% function [a,faceA,faceB]=findAB(G,interpFace,c,Kn)
% x1=G.cells.centroids(c,:)';
% theFaces=G.cells.faces(G.cells.facePos(c):G.cells.facePos(c+1)-1,1);
% myBases=interpFace.coords(theFaces,:);
% myBases=bsxfun(@minus,myBases,x1');
% myNorm=sqrt(dot(myBases,myBases,2));
% myBases=bsxfun(@rdivide,myBases,myNorm);
% Kn_norm=norm(Kn);
% Kn_unit=Kn/Kn_norm;
% myangles=bsxfun(@times,myBases,Kn_unit');
% myangles=sum(myangles,2);
% myangles=acos(myangles);
% [~,I]=sort(myangles);
% theFaces=theFaces(I);
% myBases=myBases(I,:);
% myNorm=myNorm(I);
% nf=numel(theFaces);
% flag=0;
% 
% myIndex=zeros(nf*(nf-1)/2,2);
% myCoeff=myIndex;counter=1;
% for i=1:nf-1
%     tA=myBases(i,:)';
%     tA_norm=myNorm(i);
%     for j=i+1:nf
%         tB=myBases(j,:)';
%         tB_norm=myNorm(j);
%         if(abs(det([tA tB]))>1e-9)
%             temp_a=[tA tB]\(Kn_unit);
%             temp_a(abs(temp_a)<1e-9)=0;
%             if(all(temp_a>=0))
%                 if(all(temp_a<=1))
%                     faceA=theFaces(i);
%                     faceB=theFaces(j);
%                     a=temp_a;
%                     a(1)=a(1)*Kn_norm/tA_norm;
%                     a(2)=a(2)*Kn_norm/tB_norm;
%                     flag=1;break;
%                 else
%                     myIndex(counter,:)=[i,j];
%                     myCoeff(counter,:)=temp_a;
%                     counter=counter+1;
%                 end
%             end
%         end
%     end
%     if(flag),break;end
% end
% if(~flag&&counter>1)
%     myIndex(counter:end,:)=[];myCoeff(counter:end,:)=[];
%     maxCoeff=max(myCoeff,[],2);
%     [~,ind]=min(maxCoeff);
%     i=myIndex(ind,1);j=myIndex(ind,2);
%     a=myCoeff(ind,:);
%     faceA=theFaces(i);faceB=theFaces(j);
%     tA_norm=myNorm(i);
%     tB_norm=myNorm(j);
%     a(1)=a(1)*Kn_norm/tA_norm;
%     a(2)=a(2)*Kn_norm/tB_norm;
%     flag=1;
% elseif(~flag&&counter==1)
%     for i=1:nf-1
%         tA=myBases(i,:)';
%         tA_norm=myNorm(i);
%         for j=i+1:nf
%             tB=myBases(j,:)';
%             tB_norm=myNorm(j);
%             if(abs(det([tA tB]))>1e-9)
%                 temp_a=[tA tB]\(Kn_unit);
%                 temp_a(abs(temp_a)<1e-9)=0;
%                 if(sum(temp_a)>=0)
%                     faceA=theFaces(i);
%                     faceB=theFaces(j);
%                     a=temp_a;
%                     a(1)=a(1)*Kn_norm/tA_norm;
%                     a(2)=a(2)*Kn_norm/tB_norm;
%                     flag=1;break;
%                 end
%             end
%         end
%         if(flag),break;end
%     end
% end
% 
% if(~flag)
%     error(['Decomposition failed for cell ',num2str(c)]);
% end
% end
% 
% function [a,faceA,faceB,faceC]=findABC(G,interpFace,c,Kn)
% x1=G.cells.centroids(c,:)';
% theFaces=G.cells.faces(G.cells.facePos(c):G.cells.facePos(c+1)-1,1);
% myBases=interpFace.coords(theFaces,:);
% myBases=bsxfun(@minus,myBases,x1');
% myNorm=sqrt(dot(myBases,myBases,2));
% myBases=bsxfun(@rdivide,myBases,myNorm);
% Kn_norm=norm(Kn);
% Kn_unit=Kn/Kn_norm;
% myangles=bsxfun(@times,myBases,Kn_unit');
% myangles=sum(myangles,2);
% myangles=acos(myangles);
% [~,I]=sort(myangles);
% theFaces=theFaces(I);
% myBases=myBases(I,:);
% myNorm=myNorm(I);
% nf=numel(theFaces);
% flag=0;
% 
% myIndex=zeros(nf*(nf-1)*(nf-2)/6,3);
% myCoeff=myIndex;counter=1;
% for i=1:nf-2
%     tA=myBases(i,:)';
%     tA_norm=myNorm(i);
%     for j=i+1:nf-1
%         tB=myBases(j,:)';
%         tB_norm=myNorm(j);
%         for k=j+1:nf
%             tC=myBases(k,:)';
%             tC_norm=myNorm(k);
%             if(abs(det([tA tB tC]))>1e-9)
%                 temp_a=[tA tB tC]\(Kn_unit);
%                 temp_a(abs(temp_a)<1e-9)=0;
%                 if(all(temp_a>=0))
%                     if(all(temp_a<=1))
%                         faceA=theFaces(i);
%                         faceB=theFaces(j);
%                         faceC=theFaces(k);
%                         a=temp_a;
%                         a(1)=a(1)*Kn_norm/tA_norm;
%                         a(2)=a(2)*Kn_norm/tB_norm;
%                         a(3)=a(3)*Kn_norm/tC_norm;
%                         flag=1;break;
%                     else
%                         myIndex(counter,:)=[i,j,k];
%                         myCoeff(counter,:)=temp_a;
%                         counter=counter+1;
%                     end
%                 end
%             end
%         end
%         if(flag),break;end
%     end
%     if(flag),break;end
% end
% if(~flag&&counter>1)
%     myIndex(counter:end,:)=[];myCoeff(counter:end,:)=[];
%     maxCoeff=max(myCoeff,[],2);
%     [~,ind]=min(maxCoeff);
%     i=myIndex(ind,1);j=myIndex(ind,2);k=myIndex(ind,3);
%     a=myCoeff(ind,:);
%     faceA=theFaces(i);faceB=theFaces(j);faceC=theFaces(k);
%     tA_norm=myNorm(i);
%     tB_norm=myNorm(j);
%     tC_norm=myNorm(k);
%     a(1)=a(1)*Kn_norm/tA_norm;
%     a(2)=a(2)*Kn_norm/tB_norm;
%     a(3)=a(3)*Kn_norm/tC_norm;
%     flag=1;
% elseif(~flag&&counter==1)
%     for i=1:nf-2
%         tA=myBases(i,:)';
%         tA_norm=myNorm(i);
%         for j=i+1:nf-1
%             tB=myBases(j,:)';
%             tB_norm=myNorm(j);
%             for k=j+1:nf
%                 tC=myBases(k,:)';
%                 tC_norm=myNorm(k);
%                 if(abs(det([tA tB tC]))>1e-9)
%                     temp_a=[tA tB tC]\(Kn_unit);
%                     temp_a(abs(temp_a)<1e-9)=0;
%                     if(sum(temp_a)>0)
%                         faceA=theFaces(i);
%                         faceB=theFaces(j);
%                         faceC=theFaces(k);
%                         a=temp_a;
%                         a(1)=a(1)*Kn_norm/tA_norm;
%                         a(2)=a(2)*Kn_norm/tB_norm;
%                         a(3)=a(3)*Kn_norm/tC_norm;
%                         flag=1;break;
%                     end
%                 end
%             end
%             if(flag),break;end
%         end
%         if(flag),break;end
%     end
% end
% if(~flag)
%     error(['Decomposition failed for cell ',num2str(c)]);
% end
% end


