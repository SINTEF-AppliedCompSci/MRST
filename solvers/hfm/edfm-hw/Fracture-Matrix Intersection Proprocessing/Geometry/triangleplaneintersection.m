function [truth,points,type]=triangleplaneintersection(T1,normal,point,tol)
% This function determines if a triangle and plane intersect. If so, the 
% intersection between 3D triangle T1 and the plane are provided.
% T1 is specified as a 3x3 matrix with each row containing the x,y,z 
% coordinates of T1's vertices. normal and point are used to describe the
% plane.
%
% The line output will be a matrix with 0-3 rows containing x,y,z
% coordinates of the endpoints of the intersection. Truth will be a
% true or false statement. Type will be 'triangle' 'line', 'point' or 
% 'none', depending on intersection type.



%% Determine which side T1's vertices are respective to the plane normal. 
% The possible cases are:
% 1. 2 vertices on positive side, 1 on negative
% 2. 2 vertices on positive side, 1 on plane
% 3. 1 vertex on positive, 1 on plane, 1 on negative
% 4. 1 vertex on positive side, 2 on negative
% 5. 1 vertex on plane, 2 on negative
% 6. 1 vertex on plane, 2 on positive
% 7. 2 vertices on plane, 1 on negative
% 8. 2 vertices on plane, 1 on positive
% 9. 3 vertices on positive side
% 10. 3 vertices on negative side
% 11. 3 vertices on plane
side=dot(T1-repmat(point,3,1),repmat(normal,3,1),2);
positiveside=T1(side>tol,:);
negativeside=T1(side<-tol,:);
onplane=T1(abs(side)<tol,:);
poseqside=[positiveside;onplane]; % lump positive and on plane vertices together

% Case 11
if size(onplane,1)==3
    points=onplane;
    truth=true;
    type='triangle';
    return;
end

% Cases 9 and 10
if size(positiveside,1)==3 || size(negativeside,1)==3
    points=[];
    truth=false;
    type='none';
    return;
end

% Cases 7 and 8
if size(onplane,1)==2
    points=onplane;
    truth=true;
    type='line';
    return;
end

% cases 5 and 6
if size(onplane,1)==1 && (size(positiveside,1)==2 || size(negativeside,1)==2)
    points=onplane;
    truth=true;
    type='point';
    return;
end


%% Calculate intersection line of T1 and plane for all other cases
% Assume two points on each side. If there is only one point on a side, the
% point will be repeated.
pside_1=poseqside(1,:);
pside_2=poseqside(end,:);

nside_1=negativeside(1,:);
nside_2=negativeside(end,:);

rI1=dot(normal,(point-nside_1))/dot(normal,(pside_1-nside_1));
Point1=nside_1+rI1*(pside_1-nside_1);

rI2=dot(normal,(point-nside_2))/dot(normal,(pside_2-nside_2));
Point2=nside_2+rI2*(pside_2-nside_2);

points=[Point1;Point2];
points=uniquetol(points,tol,'ByRows',true,'DataScale',1);
truth=true;

if size(points,1)==1
    type='point';
else
    type='line';
end

end
