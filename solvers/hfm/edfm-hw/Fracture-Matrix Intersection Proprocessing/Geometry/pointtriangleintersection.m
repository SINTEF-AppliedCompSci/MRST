function [truth,location] = pointtriangleintersection(point,trinodes,tol)
% POINT TRIANGLE INTERSECTION
% Detects if a point intersects a triangle and where. Inputs are point
% (containing xyz coordinates in row form) and trinodes (containing nodal
% xyz coordinates in each row). Outputs are truth (true/false) and location
% ('perimeter', 'interior', 'none')

refnode=trinodes(1,:);
node1=trinodes(2,:);
node2=trinodes(3,:);

dir1=node1-refnode;
dir2=node2-refnode;
ndir=cross(dir1,dir2); % unnecessary to make into unit normal

relposition=point-refnode;

ndist=abs(dot(relposition,ndir));
if ndist>=tol % point is not in plane
    truth=false; location='none';
    return;
end

A=[dot(dir1,dir1),dot(dir1,dir2);dot(dir2,dir1),dot(dir2,dir2)];
B=[dot(relposition,dir1);dot(relposition,dir2)];
sol=A\B;
loc1=sol(1);loc2=sol(2);

if loc1<1 && loc1>0 && loc2<1 && loc2>0 && (loc1+loc2)<1 % interior
    truth=true; location='interior';
elseif (abs(loc1)<tol && loc2<1 && loc2>0) || (abs(loc2)<tol && loc1<1 && loc1>0)|| abs(loc1+loc2-1)<tol
    % perimeter
    truth=true; location='perimeter';
else
    truth=false; location='none';
end


end

