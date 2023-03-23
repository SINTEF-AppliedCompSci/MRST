function [truth,point,location,type]=linesegmentrectangleintersect_obs(nodes,endpts,tol)
% LINE SEGMENT VS RECTANGLE INTERSECT
% Determines the intersection between a 3D line segment and a 3D rectangle.
% Inputs are nodal coordinates of the rectangle and the endpoints of the
% line segment.
% Outputs are a truth statement (true/false), point (either one or two
% points depending on the intersection type), location ('interior', 'perimeter'
% or 'none') and type ('point', 'line', 'none')

% METHODOLOGY - use function 'linesegmentplaneintersect' to determine if
% the line segment intersects an infinite plane that passes through the
% rectangle. There can be 3 outcomes:
% 1. no intersection with infinite plane
% 2. one intersection with infinite plane
% 3. two intersections with infinite plane (line)

dist=zeros(4,1);
for i=2:4
    dist(i)=norm(nodes(i,:)-nodes(1,:));
end
[~,index]=sort(dist);
dir1=nodes(index(2),:)-nodes(1,:);
dir2=nodes(index(3),:)-nodes(1,:);

planenormal=cross(dir1,dir2);
planepoint=nodes(1,:);

[intersect,xpoint,~] = linesegmentplaneintersect(endpts,planepoint,planenormal,tol);

%% NO INTERSECTION WITH INFINITE PLANE
% if the line segment has no intersection with the plane, it certainly does
% not intersect with the rectangle. Quit function.
if ~intersect
    truth=false;
    type='none';
    location='none';
    point=[];
    return;
end

%% ONE INTERSECTION WITH INFINITE PLANE
% if the line segment has one intersection with the plane, the point could
% be outside the rectangle, on the perimeter or inside the rectangle
if size(xpoint,1)==1
   relposition=xpoint-nodes(1,:);
   loc1=dot(relposition,dir1);
   loc2=dot(relposition,dir2);
   if loc1<1 && loc1>0 && loc2<1 && loc2>0 % point is within rectangle
       truth=true;
       type='point';
       location='interior';
       point=xpoint;
   elseif (abs(loc1-1)<tol || abs(loc1)<tol) && (abs(loc2-1)<tol || abs(loc2)<tol)
       truth=true;
       type='point';
       location='perimeter';
       point=xpoint;
   else
       truth=false;
       type='none';
       location='none';
       point=[];
   end 
   return;
end

%% TWO INTERSECTIONS WITH INFINITE PLANE
% if the line segment has two intersections with the plane (line segment
% is in plane), we could have one of the following scenarios:
% 1. line does not intersect rectangle
% 2. intersection line is inside the rectangle
% 3. intersection line is on the edge (perimeter)
% 4. only point intersection on the edge/vertex (perimeter)
if size(xpoint,1)==2
    point=[];
    
    %Append endpoints that are inside rectangle to point set.
    for i=1:2
        relposition=xpoint(i,:)-nodes(1,:);
        loc1=dot(relposition,dir1);
        loc2=dot(relposition,dir2);
        if loc1<1 && loc1>0 && loc2<1 && loc2>0
            point=[point;xpoint(i,:)];
        end        
    end 
    
    %Find all possible intersections with edges of rectangle, append to
    %point set.
    endpts1=xpoint;
    
    %edge1
    endpts2=nodes(index([1 2]),:);
    [~,newpoint1] = lineseglinesegintersect(endpts1,endpts2,tol);
    point=[point;newpoint1];
    
    %edge2
    endpts2=nodes(index([1 3]),:);
    [~,newpoint2] = lineseglinesegintersect(endpts1,endpts2,tol);    
    point=[point;newpoint2];
    
    %edge3
    endpts2=nodes(index([2 4]),:);
    [~,newpoint3] = lineseglinesegintersect(endpts1,endpts2,tol);  
    point=[point;newpoint3];
    
    %edge4
    endpts2=nodes(index([3 4]),:);
    [~,newpoint4] = lineseglinesegintersect(endpts1,endpts2,tol);  
    point=[point;newpoint4];
    
    %Remove repetitions in the point set
    point=uniquetol(point,tol,'ByRows',true);
    
    %Determine which of the four cases the intersection falls under
    if size(point,1)==0 % case 1: no point of intersection with rectangle
        truth=false;
        location='none';
        type='none';
    elseif size(point,1)==1 % case 4: one point of intersection on perimeter
        truth=true;
        location='perimeter';
        type='point';
    elseif size(newpoint1,1)==2 || size(newpoint2,1)==2 || size(newpoint3,1)==2 || size(newpoint4,1)==2
        % case 3: intersection line is on the perimeter
        truth=true;
        location='perimeter';
        type='line';
    else
        % case 2: intersection line is inside the rectangle
        truth=true;
        location='interior';
        type='line';
    end
return;
end

end

