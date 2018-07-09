function [truth,type] = checkplaneAABBintersect(nodes,planepoint,planenormal,tol,varargin)
% PLANE AABB INTERSECTION CHECKER - Checks if a plane and an AABB
% intersect each other. The inputs are the nodes of the AABB (row
% containing xyz coordinates), a point on the plane (row vector with xyz
% coordinates) and a normal vector to the plane (row vector with xyz
% direction coefficients). The output is a simple true or false statement.
% If requested, the type of intersection can also be returned. The possible
% types are 'vertex', 'edge', 'face' and 'interior'.
%
% Reference: Essential Mathematics for Games & Interactive Applications: A 
% Programmer's Guide by Berth & Bishop

opt=struct('positiveforface',false);
opt=merge_options(opt,varargin{:});

% Find diagonal pairs and vectors in each diagonal
[nodepair,diagdir] = diagnodes(nodes,tol);

% Find diagonal vector closest to plane normal. This is done by projecting
% the plane normal vector on each diagonal and checking which projection
% length is the longest. The longest projection length corresponds to the
% diagonal direction that is the closest to the plane normal.
projectionlength=abs(dot(repmat(planenormal,4,1),diagdir,2));
[~,index]=max(projectionlength); % find index of largest projection length
nodepair=nodepair(index,:); % use index to find/save closest diagonal pair

% Check if the diagonal intersects the plane
endpoints=nodes(nodepair,:);
[diagintersect,~,diagtype] = linesegmentplaneintersect(endpoints,planepoint,planenormal,tol);

% if diagonal does not intersect the plane, plane does not intersect
% AABB, function terminates here
% if diagonal intersects the plane and the type is 'interior', function
% also terminates here
if ~diagintersect
    truth=false;
    type='none';
    return;
elseif strcmp(diagtype,'interior')
    truth=true;
    type=diagtype;
    return;
else
    truth=true;
end

% If there is an intersection but not an interior type,
% determine the type of intersection. This is done by checking all nodes'
% distances to the plane to see how many are on the plane:
% 4 points: type is 'face'
% 2 points: type is 'edge'
% 1 point: type is 'vertex'

relposition=nodes-repmat(planepoint,8,1);
dist=abs(dot(relposition,repmat(planenormal,8,1),2));
numintersect=length(dist(dist<tol));

switch numintersect
    case 4
        type='face';
    case 2
        type='edge';
    case 1
        type='vertex';
end

% Cater for 'positiveforface', true only if AABB is on positive side of
% plane for face intersections.
if opt.positiveforface && strcmp(type,'face')
    centroid=sum(nodes)/8;
    direction=dot(centroid-planepoint,planenormal);
    truth=direction>0;
end
        
       
end

