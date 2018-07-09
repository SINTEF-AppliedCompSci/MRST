function [truth,location] = pointAABBintersection(point,AABBnodes,tol)
% POINT AABB INTERSECTION
% Detects if a point intersects an AABB and where. Inputs are point
% (containing xyz coordinates in row form) and nodes (containing nodal
% xyz coordinates in each row). Outputs are truth (true/false) and location
% ('boundary', 'interior', 'none')

xcoords=AABBnodes(:,1);
ycoords=AABBnodes(:,2);
zcoords=AABBnodes(:,3);
xmin=min(xcoords);
xmax=max(xcoords);
ymin=min(ycoords);
ymax=max(ycoords);
zmin=min(zcoords);
zmax=max(zcoords);
xpoint=point(1); ypoint=point(2); zpoint=point(3);

withinxrange=(xpoint<(xmax+tol) && xpoint>(xmin-tol));
withinyrange=(ypoint<(ymax+tol) && ypoint>(ymin-tol));
withinzrange=(zpoint<(zmax+tol) && zpoint>(zmin-tol));
onxmin=(abs(xpoint-xmin)<tol);
onxmax=(abs(xpoint-xmax)<tol);
onymin=(abs(ypoint-ymin)<tol);
onymax=(abs(ypoint-ymax)<tol);
onzmin=(abs(zpoint-zmin)<tol);
onzmax=(abs(zpoint-zmax)<tol);

if withinxrange && withinyrange && withinzrange
    truth=true;
    if onxmin || onxmax || onymin || onymax || onzmin || onzmax % if within range and on at least one of the boundaries
        location='boundary'; return;
    end
    location='interior'; return;
else
    truth=false;
    location='none';
end

end