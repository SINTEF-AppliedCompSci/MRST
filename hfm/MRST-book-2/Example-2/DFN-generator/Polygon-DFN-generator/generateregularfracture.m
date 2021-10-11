function fracture=generateregularfracture(numpoints,size,center,normal,varargin)
% this function generates a fracture that takes the shape of a regular
% polygon. The number of vertices, polygon center coordinate, normal
% direction are required as inputs. pn/pv pair 'randomize'/logic decides
% whether or not to randomize the orientation of the fracture around the
% normal direction. 'triangulate'/logic decides if triangulation is needed.
%
% the input parameter size is the distance from center to vertex. However,
% if the pn/pv pair 'circle'/logic is true, then the parameter size refers
% to radius. And numpoints will be the numver of vertices of a polygon used
% to represent the circle. In this case, the size of the regular polygon
% will be calculated from the size of a circle based on an area
% conservation principle.

opt=struct('randomize',true,'triangulate',true,'circle',false);
opt=merge_options(opt,varargin{:});

normal=normal/norm(normal); % make unit

angle=2*pi/numpoints;

if opt.circle % special case for pn/pv == 'circle'/true
%     sectorarea=size*size*angle; % removed 0.5*, compensated in next line
%     size=sqrt(sectorarea/sin(angle)); % removed 2*, compensated in previous line
    polysize=realsqrt(angle/sin(angle))*size;
else
    polysize=size;
end

rotmat=[cos(angle) -sin(angle); sin(angle) cos(angle)]; % rotation matrix

points=zeros(numpoints,2); points(1,:)=[polysize 0]; % setup points in x-y plane

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

fracture.points=points;
fracture.numpoints=numpoints;
fracture.center=center;
fracture.normal=normal;
fracture.size=size;
% fracture.circle=opt.circle;



%% triangulation
if opt.triangulate
    numtriangles=numpoints-2;
    trilist=[ones(numtriangles,1),(2:(numpoints-1))',(3:numpoints)'];
    fracture.TriangleList=trilist;   
end

end