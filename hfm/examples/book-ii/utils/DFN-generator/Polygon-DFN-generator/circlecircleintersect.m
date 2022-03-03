function truth=circlecircleintersect(circle1,circle2,tol)
% 3D circle vs circle intersection checker. Takes in circles in struct form
% with fields center, normal and radius. Tol is a
% numerical tolerance for double precision floating point numbers

% METHODOLOGY:
% 1. Check if normal directions are the same
% 2. If normal directions same, check if centers are in the same plane.
%   a. If same plane, check if distance against sum of radii.
%   b. If different plane, no intersection.
% 3. If normal directions different
%   a. Find intersection line between planes
%   b. Calculate normal distance of intersection line to circle centers
%   c. If normal distances exceed one/both circles' radii, no intersection
%   d. If normal distances within both radii, find xsect line for circles
%   e. Check if lines overlap. If no overlap, no intersection.
% 4. Pass all tests? Intersection.

normal1=circle1.normal; normal1=normal1/norm(normal1);
normal2=circle2.normal; normal2=normal2/norm(normal2);

cent1=circle1.center;
cent2=circle2.center;

rad1=circle1.radius;
rad2=circle2.radius;

%% Parallel normals

if abs(abs(dot(normal1,normal2))-1)<tol % both normals are parallel
    % Check if centers are in same plane
    normalcentdist=abs(dot(cent2-cent1,normal1)); % center distance along normal
    if normalcentdist>tol % not same plane
        truth=false; return; % no intersection, exit function
    end
    
    % exit if block if normal distance is effectively zero
    centdist=norm(cent2-cent1); % distance between center
    
    if (centdist-(rad1+rad2))>tol % center to center distance larger than sum of radii
        truth=false; return; % no intersection, exit function
    end
    
    % if successfully exit above if block, center to center distance is
    % smaller than sum of radii
    truth=true; return; % two circles intersect, exit function
end

% if does not enter above if block, the normal directions are not parallel

%% Non parallel normals
% find line of intersection between planes
[linepoint,linedir]=planeplaneintersect(cent1,normal1,cent2,normal2);

v=linedir; vdotv=dot(v,v);
% find line center distances for circle 1
w=cent1-linepoint;
closepoint1=linepoint+(dot(w,v)/vdotv)*v;
linecentdist1=norm(cent1-closepoint1);
if (linecentdist1-rad1)>tol, truth=false; return; end % no intersection if line distance is larger than radius, exit function

% find line center distances for circle 2
w=cent2-linepoint;
closepoint2=linepoint+(dot(w,v)/vdotv)*v;
linecentdist2=norm(cent2-closepoint2);
if (linecentdist2-rad2)>tol, truth=false; return; end % no intersection if line distance is larger than radius, exit function

% if did not exit function at this point, plane intersection line crosses
% both circles. We now extract line segments for both circles.
clength1=sqrt((rad1^2)-(linecentdist1^2)); % half chord length for circle 1
clength2=sqrt((rad2^2)-(linecentdist2^2)); % half chord length for circle 2

lineseg1=[(closepoint1-clength1*v);(closepoint1+clength1*v)];
lineseg2=[(closepoint2-clength2*v);(closepoint2+clength2*v)];

% we overlap the linesegments to find their intersection
relpos1=dot(lineseg1-repmat(linepoint,2,1),repmat(v,2,1),2);
relpos2=dot(lineseg2-repmat(linepoint,2,1),repmat(v,2,1),2);

overlapdist=max(0,min(max(relpos1),max(relpos2))-max(min(relpos1),min(relpos2)));

if overlapdist>tol
    truth=true; return;
end

truth=false; return;

end