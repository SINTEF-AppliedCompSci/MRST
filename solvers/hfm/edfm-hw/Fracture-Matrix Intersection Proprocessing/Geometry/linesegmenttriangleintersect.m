function [truth,point,location,type] = linesegmenttriangleintersect(nodes,endpts,tol)
% LINE SEGMENT VS TRIANGLE INTERSECT
% Determines the intersection between a 3D line segment and a 3D triangle.
% Inputs are nodal coordinates of the triangle and the endpoints of the
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

dir1=nodes(2,:)-nodes(1,:);
dir2=nodes(3,:)-nodes(1,:);

planenormal=cross(dir1,dir2);
planepoint=nodes(1,:);

[intersect,xpoint,~] = linesegmentplaneintersect(endpts,planepoint,planenormal,tol);

%% NO INTERSECTION WITH INFINITE PLANE
% if the line segment has no intersection with the plane, it certainly does
% not intersect with the triangle. Quit function.
if ~intersect
    truth=false;
    type='none';
    location='none';
    point=[];
    return;
end

%% ONE INTERSECTION WITH INFINITE PLANE
% if the line segment has one intersection with the plane, the point could
% be outside the triangle, on the perimeter or inside the triangle
if size(xpoint,1)==1
   relposition=xpoint-nodes(1,:);
   A=[dot(dir1,dir1),dot(dir1,dir2);dot(dir2,dir1),dot(dir2,dir2)];
   B=[dot(relposition,dir1);dot(relposition,dir2)];
   sol=A\B;
   loc1=sol(1);loc2=sol(2);
   if loc1<1 && loc1>0 && loc2<1 && loc2>0 && (loc1+loc2)<1 % point is within triangle
       truth=true;
       type='point';
       location='interior';
       point=xpoint;
   elseif (loc1<=1 && loc1>=0 && abs(loc2)<tol) || (loc2<=1 && loc2>=0 && abs(loc1)<tol) || abs(loc1+loc2-1)<tol
       % point is on perimeter
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
% 1. line does not intersect triangle
% 2. intersection line is inside the triangle
% 3. intersection line is on the edge (perimeter)
% 4. only point intersection on the edge/vertex (perimeter)
if size(xpoint,1)==2
    point=-1*ones(5,3); % generate 20 points to cater for excess
    count=0;
    
    %Append endpoints that are inside rectangle to point set.
    for i=1:2
        [xsect,location]=pointtriangleintersection(xpoint(i,:),nodes,tol);
        if xsect && strcmp(location,'interior')
            count=count+1;
            point(count,:)=xpoint(i,:);
        end
    end 
    
    %Find all possible intersections with edges of triangle, append to
    %point set.
    endpts1=xpoint;
    
    %edge1
    endpts2=nodes(1:2,:);
    [~,newpoint1] = lineseglinesegintersect(endpts1,endpts2,tol);
    numnewpoint=size(newpoint1,1);
    if numnewpoint>0
        point((count+1):(count+numnewpoint),:)=newpoint1;
        count=count+numnewpoint;
    end
    
    %edge2
    endpts2=nodes([1 3],:);
    [~,newpoint2] = lineseglinesegintersect(endpts1,endpts2,tol);    
    numnewpoint=size(newpoint2,1);
    if numnewpoint>0
        point((count+1):(count+numnewpoint),:)=newpoint2;
        count=count+numnewpoint;
    end
    
    %edge3
    endpts2=nodes([2 3],:);
    [~,newpoint3] = lineseglinesegintersect(endpts1,endpts2,tol);  
    numnewpoint=size(newpoint3,1);
    if numnewpoint>0
        point((count+1):(count+numnewpoint),:)=newpoint3;
        count=count+numnewpoint;
    end
        
    %Remove repetitions in the point set
    point=removeexcess(point,-1);
    point=uniquetol(point,tol,'ByRows',true,'DataScale',1);
    
    %Determine which of the four cases the intersection falls under
    if size(point,1)==0 % case 1: no point of intersection with rectangle
        truth=false;
        location='none';
        type='none';
    elseif size(point,1)==1 % case 4: one point of intersection on perimeter
        truth=true;
        location='perimeter';
        type='point';
    elseif size(newpoint1,1)==2 || size(newpoint2,1)==2 || size(newpoint3,1)==2
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

