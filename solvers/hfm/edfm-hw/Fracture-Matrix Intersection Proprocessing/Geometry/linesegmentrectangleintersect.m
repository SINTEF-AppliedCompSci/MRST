function [truth,points,location,type] = linesegmentrectangleintersect(nodes,endpts,tol)
% LINE SEGMENT VS RECTANGLE INTERSECT
% Determines the intersection between a 3D line segment and a 3D rectangle.
% Inputs are nodal coordinates of the rectangle and the endpoints of the
% line segment.
% Outputs are a truth statement (true/false), points (either one or two
% points depending on the intersection type), location ('interior', 'perimeter'
% or 'none') and type ('point', 'line', 'none')

% METHODOLOGY: We split the rectangle into two triangular segments and
% exploit the function linesegmenttriangleintersect. The possible type
% combinations for both segments are as follows:
% 1. none vs none (no intersection)
% 2. none vs point (point intersection)
% 3. none vs line (line intersection)
% 3. point vs point (point intersection)
% 4. point vs line (line intersection) - special treatment needed
% 5. line vs line (line intersection) - special treatment needed
%
% To find the location of the intersection. We use lineseglinesegintersect.
% We take all edges of the rectangle and perform intersections with the
% intersection points from the previous section. If all points are on any
% one edge, the location is on the perimeter. If no edges contain all
% points, then the location is on the interior

% figure out the edges of the rectangle by obtaining all node-node
% distances and eliminating the 2 pairs with largest distances.
possibleedges=[1 2; 1 3; 1 4; 2 3; 2 4; 3 4]; 
dist=zeros(6,1);
for i=1:6
    testpts=nodes(possibleedges(i,:),:);
    dist(i)=norm(testpts(1,:)-testpts(2,:));
end
[~,index]=sort(dist);
edges=possibleedges(index(1:4),:); % use smallest 4 distances to identify edges
T1edges=possibleedges(index([1 3]),:);
T2edges=possibleedges(index([2 4]),:);

T1=nodes(unique(T1edges(:)),:); % Segment1
T2=nodes(unique(T2edges(:)),:); % Segment2

[~,points1,~,type1] = linesegmenttriangleintersect(T1,endpts,tol);
[~,points2,~,type2] = linesegmenttriangleintersect(T2,endpts,tol);

%% Determining intersection type and intersection points
points=uniquetol([points1;points2],tol,'ByRows',true,'DataScale',1);

% special treatment for line-line combination
% point-line combination is also subject to this because using uniquetol,
% sometimes we still get 3 points.
if (strcmp(type1,'line') && strcmp(type2,'line') && size(points,1)==3) || ...
        (strcmp(type1,'point') && strcmp(type2,'line') && size(points,1)==3) ||  ...
     (strcmp(type1,'line') && strcmp(type2,'point') && size(points,1)==3)
    % after combining points, we will have a line defined by 3 points. We
    % want to remove the middle point. We do this by picking the two points
    % that give the maximum line length.
    combinations=[1 2;2 3;1 3];
    dist=zeros(3,1);
    for i=1:3
        testpts=points(combinations(i,:),:);
        dist(i)=norm(testpts(1,:)-testpts(2,:));
    end
    [~,index]=max(dist);
    maxcombi=combinations(index,:);
    points=points(maxcombi,:);  
end

if size(points,1)==0
    truth=false; type='none'; location='none'; return;
elseif size(points,1)==1
    truth=true; type='point';
elseif size(points,1)==2
    truth=true; type='line';
else
%     type1
%     type2
%     points
%     points1
%     points2
    error('Not getting a physically correct number of intersection points. Check code.');
end

%% Determining intersection location

% intersection is a point
onedge=zeros(4,1);
for i=1:4
    edge=edges(i,:);
    edgepts=nodes(edge,:);
    
    if size(points,1)==1 % intersection is a point
        edgelength=norm(edgepts(1,:)-edgepts(2,:));
        dist1=norm(points-edgepts(1,:));
        dist2=norm(points-edgepts(2,:));
        
        if abs(dist1+dist2-edgelength)<tol % point is on edge
            onedge(i)=1;
        end
    end
    
    if size(points,1)==2 % intersection is line
        [~,xpoints] = lineseglinesegintersect(edgepts,points,tol);
        
        if size(xpoints,1)==2 % points are all on edge
            onedge(i)=1;
        end
    end
    
end

if any(onedge==1)
    location='perimeter';
else
    location='interior';
end


end

