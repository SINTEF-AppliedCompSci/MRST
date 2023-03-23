function [truth,points,type] = trianglerectangleintersect(Tnodes,Rnodes,tol)
% This function determines if a triangle and rectangle intersect. It also 
% returns, if requested, the intersection of the two 3D planes
% The inputs are Tnodes and Rnodes. These are 3x3 and 4x3
% matrices with each row being the x,y,z coordinates of the vertices of the
% objects. The outputs are truth (true/false), points (matrix with each
% row containing vector position of intersection point), type ('polygon',
% 'line', 'point', 'none')

% METHODOLOGY: Adapt from triangletriangleintersection. Compare triangle
% with plane. Then rectangle with plane. Then perform an intersection of
% the tp and rp intersections.

%% Define infinite planes
% Plane for triangle
V1_T=Tnodes(1,:);
V2_T=Tnodes(2,:);
V3_T=Tnodes(3,:);
normal_T=cross(V2_T-V1_T,V3_T-V1_T)/norm(cross(V2_T-V1_T,V3_T-V1_T));

% Plane for rectangle
V1_R=Rnodes(1,:);
V2_R=Rnodes(2,:);
V3_R=Rnodes(3,:);
normal_R=cross(V2_R-V1_R,V3_R-V1_R)/norm(cross(V2_R-V1_R,V3_R-V1_R));

%% Check triangle vs plane and rectangle vs plane intersections
[TxRP,pointsTxRP,typeTxRP]=triangleplaneintersection(Tnodes,normal_R,V1_R,tol);
[RxTP,pointsRxTP,typeRxTP]=rectangleplaneintersection(Rnodes,normal_T,V1_T,tol);

%% Combine the results of the above triangle/rectangle-plane intersections
% The output from the above can be the following types of combinations
% 1. empty set with any other type (no intersection at all)
% 2. point and point (same points mean intersect at point)
% 3. point and line segment (possible point intersection)
% 4. point and triangle/rectangle (possible point intersection)????? unrealistic
% 5. line and triangle (use linesegmenttriangleintersect)?????? unrealistic
% 6. triangle and rectangle (both in same plane)
% 7. line and line (use lineseglinesegintersect)
% 8. line and rectangle (use linesegmentrectangleintersect)?????? unrealistic

% case 1: empty set vs all
if ~TxRP || ~RxTP
    truth=false;points=[];type='none';
    return;
end

% case 2: point vs point
if strcmp(typeTxRP,'point') && strcmp(typeRxTP,'point')
   points=uniquetol([pointsTxRP;pointsRxTP],tol,'ByRows',true,'DataScale',1);
   if size(points,1)==2 % the two points don't coincide
       truth=false;points=[];type='none';
   elseif size(points,1)==1 % the two points coincide
       truth=true;type='point';
   end
   return;
end

% case 3: point vs line
if (strcmp(typeTxRP,'point') && strcmp(typeRxTP,'line')) || (strcmp(typeTxRP,'line') && strcmp(typeRxTP,'point'))
    if strcmp(typeTxRP,'point')
        TPpoint=pointsTxRP;
        TPline=pointsRxTP;
    else
        TPpoint=pointsRxTP;
        TPline=pointsTxRP;
    end
    
    linelength=norm(TPline(1,:)-TPline(2,:));
    pointendpt1dist=norm(TPpoint-TPline(1,:));
    pointendpt2dist=norm(TPpoint-TPline(2,:));
    
    % point is on line if the following parameter is 0
    difflength=pointendpt1dist+pointendpt2dist-linelength;
    
    if abs(difflength)<tol
        truth=true;points=TPpoint;type='point';
    else
        truth=false;points=[];type='none';
    end
    return;
end

% case 7: line and line
if strcmp(typeTxRP,'line') && strcmp(typeRxTP,'line')
    line1=pointsTxRP;
    line2=pointsRxTP;
       
    [truth,points] = lineseglinesegintersect(line1,line2,tol);
    
    if size(points,1)==0
        type='none';
    elseif size(points,1)==1
        type='point';
    else
        type='line';
    end
    
    return;
end

% case 6: triangle & rectangle
% Method
%   a. determine/save if rect endpoints are in tri
%   b. determine/save if tri endpoints are in rect
%   c. determine/save rect edges vs tri intersections
%   d. determine/save tri edges vs rect intersections
%   e. eliminate repetitions
%   f. determine if intersection is point, line, or polygon
if strcmp(typeTxRP,'triangle') && strcmp(typeRxTP,'rectangle')
    tri=pointsTxRP;
    rect=pointsRxTP;
    
    points=-1*ones(20,3); % preallocate 20 rows, remove later
    count=0;
    
    % step a: determine/save if rect endpoints are in tri
    for i=1:4
        [rint,location] = pointtriangleintersection(rect(i,:),tri,tol);
        if rint && strcmp(location,'interior')
            count=count+1;
            points(count,:)=rect(i,:);
        end
    end
    
    % step b: determine/save if tri endpoints are in rect
    for i=1:3
        [tinr,location] = pointrectangleintersection(tri(i,:),rect,tol);
        if tinr && strcmp(location,'interior')
            count=count+1;
            points(count,:)=tri(i,:);
        end
    end
    
    % step c: determine/save rect edges vs tri intersections
    % determine edges
    possibleedges=[1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
    dist=zeros(6,1);
    for i=1:6
        testpts=rect(possibleedges(i,:),:);
        dist(i)=norm(testpts(1,:)-testpts(2,:));
    end
    [~,index]=sort(dist);
    edges=possibleedges(index(1:4),:); % use smallest 4 distances to identify edge 
    for i=1:4
        [~,xpoints,~,~] = linesegmenttriangleintersect(tri,rect(edges(i,:),:),tol);
        numxpoints=size(xpoints,1);
        if numxpoints>0
            points((count+1):(count+numxpoints),:)=xpoints;
            count=count+numxpoints;
        end
    end
    
    % step d: determine/save tri edges vs rect intersections
    edges=[1 2; 2 3; 1 3]; % row indices of nodes that form triangle edges
    for i=1:3
        [~,xpoints,~,~] = linesegmentrectangleintersect(rect,tri(edges(i,:),:),tol);
        numxpoints=size(xpoints,1);
        if numxpoints>0
            points((count+1):(count+numxpoints),:)=xpoints;
            count=count+numxpoints;
        end
    end
    
    % step e: eliminate repetitions
    points=removeexcess(points,-1);
    points=uniquetol(points,tol,'ByRows',true,'DataScale',1);
    
    % step f: determine if intersection is point, line, or polygon
    if size(points,1)==0
        truth=false;type='none';
    elseif size(points,1)==1
        truth=true;type='point';
    elseif size(points,1)==2
        truth=true;type='line';
    else
        truth=true;type='polygon';
    end
    
    return;
end

truth=false;type='none';points=[];
% justification for last line. Sometimes we get triangular gridcells that
% are actually just a line. Not sure why this is so, but we discard these
% occurences.
% Also, we treat polygon vs non polygon occurences as numerical errors.
end

