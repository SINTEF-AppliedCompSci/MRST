function [truth,points,type]=triangletriangleintersection(T1,T2,tol)
% This function determines if triangles intersect. It also returns, if 
% requested, the intersection of two 3D triangular
% planes. The inputs are T1 and T2. These are 3x3
% matrices with each row being the x,y,z coordinates of the vertices of the
% triangles. The outputs are truth (true/false), points (matrix with each
% row containing vector position of intersection point), type ('polygon',
% 'line', 'point', 'none')

%% Check triangle vs infinite plane intersections
% Triangle 1
V1_1=T1(1,:);
V2_1=T1(2,:);
V3_1=T1(3,:);
normal_1=cross(V2_1-V1_1,V3_1-V1_1)/norm(cross(V2_1-V1_1,V3_1-V1_1));

% Triangle 2
V1_2=T2(1,:);
V2_2=T2(2,:);
V3_2=T2(3,:);
normal_2=cross(V2_2-V1_2,V3_2-V1_2)/norm(cross(V2_2-V1_2,V3_2-V1_2));

[T1xP2,pointsT1xP2,typeT1xP2]=triangleplaneintersection(T1,normal_2,V1_2,tol);
[T2xP1,pointsT2xP1,typeT2xP1]=triangleplaneintersection(T2,normal_1,V1_1,tol);

%% Combine the results of the above triangle-plane intersections
% The output from the above can be the following types of combinations
% 1. empty set with any other type (no intersection at all)
% 2. point and point (same points mean intersect at point)
% 3. point and line segment (possible point intersection)
% 4. point and triangle (possible point intersection)?????????? unrealistic
% 5. line and triangle (use linesegmenttriangleintersect)????????? unrealistic
% 6. triangle and triangle (both in same plane)
% 7. line and line (use lineseglinesegintersect)

% case 1: one empty set
if ~T1xP2 || ~T2xP1
    truth=false;points=[];type='none';
    return;
end

% case 2: point-point
if strcmp(typeT1xP2,'point') && strcmp(typeT2xP1,'point')
   points=uniquetol([pointsT1xP2;pointsT2xP1],tol,'ByRows',true,'DataScale',1);
   if size(points,1)==2 % the two points don't coincide
       truth=false;points=[];type='none';
   elseif size(points,1)==1 % the two points coincide
       truth=true;type='point';
   end
   return;
end

% case 3: point-line
if (strcmp(typeT1xP2,'point') && strcmp(typeT2xP1,'line')) || (strcmp(typeT1xP2,'line') && strcmp(typeT2xP1,'point'))
    if strcmp(typeT1xP2,'point')
        TPpoint=pointsT1xP2;
        TPline=pointsT2xP1;
    else
        TPpoint=pointsT2xP1;
        TPline=pointsT1xP2;
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

    % % case 4: point-triangle
    % if (strcmp(typeT1xP2,'point') && strcmp(typeT2xP1,'triangle')) || (strcmp(typeT1xP2,'triangle') && strcmp(typeT2xP1,'point'))
    %     if strcmp(typeT1xP2,'point')
    %         TPpoint=pointsT1xP2;
    %         TPnodes=pointsT2xP1;
    %     else
    %         TPpoint=pointsT2xP1;
    %         TPnodes=pointsT1xP2;
    %     end
    %     
    %     % determine if point and triangle intersect
    %     [PTintersect,~] = pointtriangleintersection(TPpoint,TPnodes,tol);
    %     
    %     if PTintersect
    %         truth=true;points=TPpoint;type='point';
    %     else
    %         truth=false;points=[];type='none';
    %     end
    %     return;
    % end

    % % case 5: line-triangle
    % if (strcmp(typeT1xP2,'line') && strcmp(typeT2xP1,'triangle')) || (strcmp(typeT1xP2,'triangle') && strcmp(typeT2xP1,'line'))
    %     if strcmp(typeT1xP2,'line')
    %         TPline=pointsT1xP2;
    %         TPnodes=pointsT2xP1;
    %     else
    %         TPline=pointsT2xP1;
    %         TPnodes=pointsT1xP2;
    %     end
    %     
    %     % determine if line and triangle intersect
    %     [truth,points,~,type] = linesegmenttriangleintersect(TPnodes,TPline,tol);
    %     return;
    % end

% case 6: triangle-triangle
% Method
%   a. determine/save if tri2 endpoints are in tri1
%   b. determine/save if tri1 endpoints are in tri2
%   c. determine/save tri1 edges vs tri2 intersections
%   d. determine/save tri2 edges vs tri1 intersections
%   e. eliminate repetitions
%   f. determine if intersection is point, line, or polygon
if strcmp(typeT1xP2,'triangle') && strcmp(typeT2xP1,'triangle')
    tri1=pointsT1xP2;
    tri2=pointsT2xP1;
    
    points=-1*ones(20,3); % preallocate 20 rows, remove excess later
    count=0;
    
    % step a: determine/save if tri2 endpoints are in tri1
    for i=1:3
        [t2int1,location] = pointtriangleintersection(tri2(i,:),tri1,tol);
        if t2int1 && strcmp(location,'interior')
            count=count+1;
            points(count,:)=tri2(i,:);
        end
    end
    
    % step b: determine/save if tri1 endpoints are in tri2
    for i=1:3
        [t1int2,location] = pointtriangleintersection(tri1(i,:),tri2,tol);
        if t1int2 && strcmp(location,'interior')
            count=count+1;
            points(count,:)=tri1(i,:);
        end
    end
    
    % step c: determine/save tri1 edges vs tri2 intersections
    edges=[1 2; 2 3; 1 3]; % row indices of nodes that form triangle edges
    for i=1:3
        [~,xpoints,~,~] = linesegmenttriangleintersect(tri2,tri1(edges(i,:),:),tol);
        numxpoints=size(xpoints,1);
        if numxpoints>0
            points((count+(1:numxpoints)),:)=xpoints;
            count=count+numxpoints;
        end
    end
    
    % step d: determine/save tri2 edges vs tri1 intersections
    for i=1:3
        [~,xpoints,~,~] = linesegmenttriangleintersect(tri1,tri2(edges(i,:),:),tol);
        numxpoints=size(xpoints,1);
        if numxpoints>0
            points((count+(1:numxpoints)),:)=xpoints;
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

% case 7: line-line
if strcmp(typeT1xP2,'line') && strcmp(typeT2xP1,'line')
    line1=pointsT1xP2;
    line2=pointsT2xP1;
       
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

end


