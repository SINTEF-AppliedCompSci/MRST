function points=circlepoints(numpoints,center,normal,radius,varargin)
% EQUIVALENT POLYGON GENERATOR FOR 3D CIRCLE - This function takes in the
% normal direction, center and radius of a 3D circle.
%
% Then, a polygon is generated which approximates the circle. The output
% is an array of vertex coordinates for each row.
%
% The function also requires the number of polygon vertices to be
% specified. Then the function will perform iterations to ensure that the
% area of the polygon is within tolerance of the actual elliptical area.
%
% Optional 'pn'/'pv' pair 'randomize' can be used to rotate the polygon
% randomly around the normal direction. 

opt=struct('randomize',false);
opt=merge_options(opt,varargin{:});

angle=2*pi/numpoints;
sectorarea=radius*radius*angle; % removed 0.5*, compensated in next line
pointdist=sqrt(sectorarea/sin(angle)); % removed 2*, compensated in previous line

normal=normal/norm(normal); % make unit

rotmat=[cos(angle) -sin(angle); sin(angle) cos(angle)]; % rotation matrix

points=zeros(numpoints,2); points(1,:)=[pointdist 0]; % setup points in x-y plane

for i=2:numpoints
    points(i,:)=(rotmat*points(i-1,:)')';
end

if opt.randomize
    randangle=angle*rand; % random angle between 0 and 'angle'. This choice is based on symmetry
    randrot=[cos(randangle) -sin(randangle); sin(randangle) cos(randangle)]; % random rotation matrix
    points=(randrot*points')'; % rotate points based on random angle
end

points=[points,zeros(numpoints,1)]; % add in the z-coordinate

points=rotatePlane(points,normal,'useNormal',[0 0 1]); % rotate plane to circle plane

points=points+repmat(center,numpoints,1); % translate points using circle center

end

