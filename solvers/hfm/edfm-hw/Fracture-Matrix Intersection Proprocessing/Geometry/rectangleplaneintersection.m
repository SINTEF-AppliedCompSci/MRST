function [truth,points,type]=rectangleplaneintersection(Rnodes,planenormal,planepoint,tol)
% This function determines if a rectangle and plane intersect. If so, the 
% intersection between 3D rectangle and the plane are provided.
% Rnodes is specified as a 4x3 matrix with each row containing the x,y,z 
% coordinates of R's vertices. normal and point are used to describe the
% plane.
%
% The line output will be a matrix with 0-2 or 4 rows containing x,y,z
% coordinates of the endpoints of the intersection. Truth will be a
% true or false statement. Type will be 'rectangle' 'line', 'point' or 
% 'none', depending on intersection type.

% METHODOLOGY: split rectangle into two triangles and use the
% triangleplaneintersection function to determine segment intersections.
% Then the intersections are combined together. The possible outcome
% combinations from the two segments are
% 1. none and none (no intersection)
% 2. none and point (point intersection)
% 3. none and line (line intersection)
% 4. point and point (point intersection)
% 5. point and line (line intersection) - special treatment needed
% 6. line and line (line intersection) - special treatment needed
% 7. triangle and triangle (rectangle intersection)
%
% impossible scenarios are:
% 8. none and triangle
% 9. point and triangle
% 10. line and triangle

dist=zeros(4,1);
for i=2:4
    dist(i)=norm(Rnodes(i,:)-Rnodes(1,:));
end

[~,index]=sort(dist);

T1=Rnodes(index([1 2 4]),:); % Segment1
T2=Rnodes(index([1 3 4]),:); % Segment2

[~,points1,type1]=triangleplaneintersection(T1,planenormal,planepoint,tol);
[~,points2,type2]=triangleplaneintersection(T2,planenormal,planepoint,tol);

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
    truth=false; type='none';
elseif size(points,1)==1
    truth=true; type='point';
elseif size(points,1)==2
    truth=true; type='line';
elseif size(points,1)==4
    truth=true; type='rectangle';
else
    disp(points);
    disp(type1);
    disp(type2);
    disp('Not getting a physically correct number of intersection points. Check code.');
    truth=false; type='none'; points=[];
end


end

