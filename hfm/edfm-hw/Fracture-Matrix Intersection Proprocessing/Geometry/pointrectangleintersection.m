function [truth,location] = pointrectangleintersection(point,rectnodes,tol)
% POINT RECTANGLE INTERSECTION
% Detects if a point intersects a rectangle and where. Inputs are point
% (containing xyz coordinates in row form) and rectnodes (containing nodal
% xyz coordinates in each row). Outputs are truth (true/false) and location
% ('perimeter', 'interior', 'none')

%% quick first check to see if point is in plane
node1=rectnodes(1,:);
node2=rectnodes(2,:);
node3=rectnodes(3,:);

dir2=node2-node1;
dir3=node3-node1;
ndir=cross(dir2,dir3); % unnecessary to make into unit normal

relposition=point-node1;

ndist=abs(dot(relposition,ndir));
if ndist>=tol % point is not in plane
    truth=false; location='none';
    return;
end

%% Check to see if point is in rectangle
node4=rectnodes(4,:);
dir4=node4-node1;

dist2=norm(dir2); dist3=norm(dir3); dist4=norm(dir4);
dir=[dir2;dir3;dir4];
[~,index]=sort([dist2,dist3,dist4]);

dirx=dir(index(1),:);
diry=dir(index(2),:);

loc1=dot(relposition,dirx)/(norm(dirx)^2);
loc2=dot(relposition,diry)/(norm(diry)^2);

if loc1<1 && loc1>0 && loc2<1 && loc2>0 % interior
    truth=true; location='interior';
elseif ((abs(loc1)<tol || abs(loc1-1)<tol) && loc2<1 && loc2>0) || ((abs(loc2)<tol || abs(loc2-1)<tol) && loc1<1 && loc1>0) || ...
        ((abs(loc1)<tol || abs(loc1-1)<tol) && (abs(loc2)<tol || abs(loc2-1)<tol))
    % perimeter
    truth=true; location='perimeter';
else
    truth=false; location='none';
end


end

