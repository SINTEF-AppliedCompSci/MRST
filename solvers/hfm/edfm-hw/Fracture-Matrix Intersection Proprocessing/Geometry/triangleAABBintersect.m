function [truth,points,type,location] = triangleAABBintersect(trinodes,AABBnodes,tol)
% INTERSECTION BETWEEN 3D TRIANGLE AND AABB. This function determines if a
% triangle intersects an AABB. The inputs trinodes and AABBnodes are nodes
% containing row vectors of xyz coordinates. The outputs are truth
% (true/false), points (row vectors of intersection polygon or []), type
% ('polygon', 'line', 'point', 'none'), location ('interior', 'boundary', 
% 'none')

% METHODOLOGY
% 1. Use checkplaneAABBintersect
%   a. plane does not cross AABB (quick filter)
%   b. plane touches vertex (perform point-triangle intersection)
%   c. plane touches edge (perform line-triangle intersection)
%   d. plane touches face (perform triangle-rectangle intersection)
%   e. plane crosses interior (perform full intersection check - Step 2)
% 2. Perform full intersection check
%   a. determine/save triangle vertices within AABB
%   b. determine/save triangle and AABB face intersections
%   c. remove duplicates

dir1=trinodes(2,:)-trinodes(1,:);
dir2=trinodes(3,:)-trinodes(1,:);
trinormal=cross(dir1,dir2); trinormal=trinormal/norm(trinormal);
[AABBxTP,AABBxTPtype] = checkplaneAABBintersect(AABBnodes,trinodes(1,:),trinormal,tol);

%% QUICK FILTER
% if plane does not intersect AABB, triangle certainly does not intersect
% AABB

if ~AABBxTP
    truth=false; points=[]; type='none'; location='none'; return;
end

%% PLANE INTERSECTS WITH VERTEX
if strcmp(AABBxTPtype,'vertex')
    truth=false;points=[];type='none';location='none'; % default output
    % perform point-triangle intersection (8 times)
    for i=1:8
        [VxT,~] = pointtriangleintersection(AABBnodes(i,:),trinodes,tol);
        if VxT % if intersection is found, there is no need to proceed further
            truth=true;points=AABBnodes(i,:);type='point';location='boundary';return;
        end
    end
    return;
end

%% PLANE INTERSECTS WITH EDGE
if strcmp(AABBxTPtype,'edge')
    % find edges
    edges=findAABBedges(AABBnodes,tol);
    truth=false;points=[];type='none';location='none'; % default output
    % perform line-triangle intersection (12 times)
    for i=1:12
        linei=AABBnodes(edges(i,:),:);
        [ExT,xpoints,~,xtype] = linesegmenttriangleintersect(trinodes,linei,tol);
        if ExT && strcmp(xtype,'point')
            truth=true;points=xpoints;type=xtype;location='boundary';
        elseif ExT && strcmp(xtype,'line') % if we find a line on an edge, no need to proceed further
            truth=true;points=xpoints;type=xtype;location='boundary';return;
        end
    end
    % if we successfully exit for loop, the intersection point is either
    % 'point' or 'none'
    return;
end

%% PLANE INTERSECTS WITH FACE
if strcmp(AABBxTPtype,'face')
    truth=false;points=[];type='none';location='none'; % default output
    
    % find faces
    faces=findAABBfaces(AABBnodes,tol);
    
    % perform triangle-rectangle intersection (6 times)
    for i=1:6
        facenodesi=AABBnodes(faces(i,:),:);
        [FxT,xpoints,xtype] = trianglerectangleintersect(trinodes,facenodesi,tol);
        if FxT && strcmp(xtype,'polygon') % if intersection is polygon, no need to proceed further
            truth=true;points=xpoints;type=xtype;location='boundary';return;
        elseif strcmp(type,'line') % if intersection is not a polygon and a line has been detected previously, go to next loop
            continue;
        elseif FxT && strcmp(xtype,'line') % if intersection is not a polygon, no line detected previously, save this line
            truth=true;points=xpoints;type=xtype;location='boundary';
        elseif FxT && strcmp(xtype,'point') % if intersection is not polygon nor line, no line detected previously, save this point
            truth=true;points=xpoints;type=xtype;location='boundary';
        end 
    end
    
    % if function successfully exits block, then type is either none,
    % point, or line
    
    return;
end

%% PLANE INTERSECTS WITH INTERIOR
% perform full intersection check

points=-1*ones(20,3); % preallocate 20 rows. Remove later.
count=0;

% find faces
faces=findAABBfaces(AABBnodes,tol);

% determine/save triangle vertices within AABB
for i=1:3 % 3 triangle vertices
    [VinAABB,~] = pointAABBintersection(trinodes(i,:),AABBnodes,tol);
    if VinAABB % if point is in or on AABB, save it
        count=count+1;
        points(count,:)=trinodes(i,:);
    end
end
   
% determine/save triangle and AABB face intersections (6 faces)
for i=1:6
    facenodesi=AABBnodes(faces(i,:),:);
    [FxT,xpoints,~] = trianglerectangleintersect(trinodes,facenodesi,tol);
    if FxT % if triangle intersects face, save points
        numxpoints=size(xpoints,1);
        points((count+1):(count+numxpoints),:)=xpoints;
        count=count+numxpoints;
    end
end

% eliminate duplicates
points=removeexcess(points,-1);
points=uniquetol(points,tol,'ByRows',true,'DataScale',1);

if size(points,1)==0
    truth=false;type='none';location='none';
elseif size(points,1)==1
    truth=true;type='point';location='boundary'; % impossible to have interior point only
elseif size(points,2)==2
    truth=true;type='line';location='boundary'; % impossible to have interior line only
else
    truth=true;type='polygon';location='interior';
end
    
end




