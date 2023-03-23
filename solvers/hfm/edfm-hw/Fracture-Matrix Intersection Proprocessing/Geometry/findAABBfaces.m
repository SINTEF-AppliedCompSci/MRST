function faces=findAABBfaces(AABBnodes,tol)
% FIND AABB FACES. Takes in the nodal positions of an AABB (each row
% containing xyz coordinates). The output is a 6x4 matrix with each row
% containing the row indices of vertices that form a face.

% METHODOLOGY: We make use of the fact that an AABB is axis-aligned. We
% first extract the max/min x,y,z coordinates. Then for each one of them,
% we find the four vertices containing them.

xcoords=AABBnodes(:,1);
ycoords=AABBnodes(:,2);
zcoords=AABBnodes(:,3);
xmin=min(xcoords);
xmax=max(xcoords);
ymin=min(ycoords);
ymax=max(ycoords);
zmin=min(zcoords);
zmax=max(zcoords);

faces=zeros(6,4); % initialize faces

% left face (xmin)
onface=(abs(xcoords-xmin*ones(8,1))<tol);
index=find(onface);
faces(1,:)=index';

% right face (xmax)
onface=(abs(xcoords-xmax*ones(8,1))<tol);
index=find(onface);
faces(2,:)=index';

% front face(ymin)
onface=(abs(ycoords-ymin*ones(8,1))<tol);
index=find(onface);
faces(3,:)=index';

% back face(ymax)
onface=(abs(ycoords-ymax*ones(8,1))<tol);
index=find(onface);
faces(4,:)=index';

% top face(zmin)
onface=(abs(zcoords-zmin*ones(8,1))<tol);
index=find(onface);
faces(5,:)=index';

% bottom face(zmax)
onface=(abs(zcoords-zmax*ones(8,1))<tol);
index=find(onface);
faces(6,:)=index';

end

