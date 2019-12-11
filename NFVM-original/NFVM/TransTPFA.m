function [ T ] = TransTPFA( G,rock,bc)
%Compute static transmissibility of faces using TPFA
%  Applicable to general 2-D and 3-D grids
%  G - Grid structure of MRST
%  rock - rock structure of MRST
%  bc - myself designed BC structure,function handle for pressure boundary
%  and outward flux for Neumann boundary
%  T - half-face transmissiblity including bc condition

%%
K=permTensor(rock,G.griddim);
K=reshape(K',G.griddim,G.griddim,[]);
T=zeros(G.faces.num,2);
N=G.faces.neighbors;
for i_face=1:G.faces.num  
    if(all(N(i_face,:)~=0))  % ------------------------------internal faces
        c1=N(i_face,1);c2=N(i_face,2);
        x1=G.cells.centroids(c1,:)';
        x2=G.cells.centroids(c2,:)';
        xm=G.faces.centroids(i_face,:)';
        d1=xm-x1;d2=xm-x2;
        face_normal=G.faces.normals(i_face,:)';
        perm1=K(:,:,c1);perm2=K(:,:,c2);
        T(i_face,1)=abs(d1'*perm1*face_normal/dot(d1,d1));
        T(i_face,2)=abs(d2'*perm2*face_normal/dot(d2,d2));       
    else  %------------------------------------ faces lying on the boundary
        c=max(N(i_face,:));
        xc=G.cells.centroids(c,:)';
        id=find(bc.face==i_face,1);
        xm=G.faces.centroids(i_face,:)';
        if(strcmpi('pressure',bc.type{id}))
            d=xm-xc;
            face_normal=G.faces.normals(i_face,:)';
            perm=K(:,:,c);
            T(i_face,1)=abs(d'*perm*face_normal/dot(d,d));
            T(i_face,2)=-T(i_face,1)*bc.value{id}(xm);
        else
            T(i_face,2)=bc.value{id}(xm)*G.faces.areas(i_face);
        end        
    end
end
end

