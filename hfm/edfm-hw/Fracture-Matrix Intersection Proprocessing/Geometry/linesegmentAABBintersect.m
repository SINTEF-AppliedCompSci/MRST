function [truth,points,type,location] = linesegmentAABBintersect(endpts,nodes,tol)
% CHECKER FOR AABB VS LINE SEGMENT INTERSECTION
%
% Inputs are the endpoints of the line segment (each row containing xyz
% coordinates of endpoint) and nodes of the AABB (each row containing xyz
% of the vertices).
%
% Outputs are a truth statement (true for intersection, false for no
% intersection),point (xyz coordinates of points defining line inside
% AABB or []),location ('boundary', 'interior', 'none') and type
% ('point', 'line', 'none').
%
% Reference: Essential Mathematics for Games & Interactive Applications: A 
% Programmer's Guide by Verth & Bishop

%% COLLISION DETECTION BETWEEN AABB AND LINE SEGMENT

refpoint=endpts(1,:); % endpoint1
unitdir=(endpts(2,:)-endpts(1,:))/norm(endpts(2,:)-endpts(1,:)); % unit vector pointing from endpoint1 to endpoint2
length=norm(endpts(2,:)-endpts(1,:)); % length of segment

% initialize maxS=0 & minT=length
maxS=0;
minT=length;

% loop through all 3 directions (x, y and z)
maxco=zeros(1,3);
minco=zeros(1,3);
for i=1:3
    v=unitdir(i);
    P=refpoint(i);
    maxco(i)=max(nodes(:,i)); % maximum coordinate in direction, save for use later
    minco(i)=min(nodes(:,i)); % minimum coordinate in direction, save for use later
    
    if v>=0
        s=(minco(i)-P)/v;
        t=(maxco(i)-P)/v;
    else
        s=(maxco(i)-P)/v;
        t=(minco(i)-P)/v;
    end
    
    if s>maxS
        maxS=s;
    end
    
    if t<minT
        minT=t;
    end
    
    if maxS>minT
        truth=false;
        points=[];
        type='none';
        location='none';
        return; % exit function at this point
    end
end

% if exit for loop successfully, then intersection exists
truth=true;

%% CALCULATION FOR LINE SEGMENT INSIDE AABB OR POINT ON AABB
% possible types of intersection line segments are as follows
% case 1: line segment is inside the AABB
% case 2: line segment is on boundary of the AABB
% case 3: point intersection on boundary of AABB

points=-1*ones(10,3); % initialize set of 20 points
count=0;
xmax=maxco(1); xmin=minco(1);
ymax=maxco(2); ymin=minco(2);
zmax=maxco(3); zmin=minco(3);

% determine if endpoints are inside AABB, if they are, append to points set.
for i=1:2
    x=endpts(i,1); y=endpts(i,2); z=endpts(i,3);
    if (x<xmax)&&(x>xmin)&&(y<ymax)&&(y>ymin)&&(z<zmax)&&(z>zmin)
        count=count+1;
        points(count,:)=[x,y,z];
    end
end

% determine intersection points with the six faces
facenodes=cell(6,1);

% left face (1)
facenodes{1}=[xmin ymin zmin; xmin ymax zmin; xmin ymax zmax; xmin ymin zmax];

% right face(2)
facenodes{2}=[xmax ymin zmin; xmax ymax zmin; xmax ymax zmax; xmax ymin zmax];

% front face (3)
facenodes{3}=[xmin ymin zmin; xmax ymin zmin; xmax ymin zmax; xmin ymin zmax];

% back face (4)
facenodes{4}=[xmin ymax zmin; xmax ymax zmin; xmax ymax zmax; xmin ymax zmax];

% bottom face (5)
facenodes{5}=[xmin ymin zmin; xmax ymin zmin; xmax ymax zmin; xmin ymax zmin];

% top face (6)
facenodes{6}=[xmin ymin zmax; xmax ymin zmax; xmax ymax zmax; xmin ymax zmax];

xonface=zeros(6,1); % number of intersections on each face (initialization step)
for i=1:6
    [~,xpoints,~,~] = linesegmentrectangleintersect(facenodes{i},endpts,tol);
    xonface(i)=size(xpoints,1);
    if xonface(i)>0
        points((count+1):(count+xonface(i)),:)=xpoints;
        count=count+xonface(i);
    end
end

points=removeexcess(points,-1);
points=uniquetol(points,tol,'ByRows',true,'DataScale',1);


if size(points,1)==1 % case 3: intersection point on boundary
    type='point';
    location='boundary';
elseif any(xonface==2) % case 2: intersection line on boundary
    type='line';
    location='boundary';
else % case 1: intersection line crosses interior of AABB
    type='line';
    location='interior';
end


end





