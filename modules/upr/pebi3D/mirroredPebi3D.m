function [G] = mirroredPebi3D(pts, bound)
% Creates 3D voronoi diagram inside a convex boundary
%
% SYNOPSIS:
%   G = mirroredPebi3D(pts, bound)
%
% PARAMETERS:
%   p         A nx3 array containing the voronoi sites. 
%   bound     A kx3 array containing the vertices of the bounding 
%             polyhedron. The bounding surface is interpreted as the convex
%             hull of bound. 
%
% RETURNS:
%   G                - Valid MRST grid definition.  
%
% EXAMPLE:
% dt = 0.15;
% [X,Y,Z] = ndgrid(dt:dt:1-dt);
% X(1:2:end) = X(1:2:end) + dt/2;
% Y(1:2:end) = Y(1:2:end) + dt/2;
% Z(1:2:end) = Z(1:2:end) + dt/2;
% pts = [X(:), Y(:),Z(:)];
% bound = [0,0,0; 1,0,0; 1,1,0; 0,1,0; 0,0,1; 1,0,1; 1,1,1;0,1,1];
% G = mirroredPebi3D(pts, bound);
% plotGrid(G);
% axis equal
%
% SEE ALSO:
%   compositePebiGrid2D, pebi, voronoi2mrstGrid3D, clipGrid.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}


% Calculate convex hull and remove all points outside domain.
dt = delaunayTriangulation(bound);
K = convexHull(dt);
inHull = ~isnan(pointLocation(dt, pts));
pts = pts(inHull,:);

% Calculate boundary normals
normals = calcNormals(K, bound);
% Find and remove boundary planes that are the same
remove = remParPlanes(bound(K(:,1),:), normals);
K = K(~remove,:);
normals = normals(~remove,:);

% mirror points around planes
m = size(pts,1);
p = zeros(m*(size(K,1)+1),3);
p(1:m,:) = pts;
for i = 1:size(K,1)
    p(m*i+1:m*(i+1),:) = mirror(pts, bound(K(i,1),:),normals(i,:));
end

% Generate voronoi grid
[V, C] = voronoin(p);

% Find auxillary cells
remove = false(numel(C),1);
for i = 1:numel(C)
    avgCell = sum(V(C{i},:),1)/size(C{i},2);
    if any(isinf(avgCell)) || isnan(pointLocation(dt, avgCell));
        remove(i) = true;
    end
end

C = C(~remove);
V   = V(2:end,:);
C   = cellfun(@(c) c-1, C,'un',false);
% Convert to mrst Grid-structure
G = voronoi2mrstGrid3D(V, C);    
end


function r = remParPlanes(x0, n)
r = false(size(x0,1),1);
for i = 1:size(n,1)-1
    % test if point is on plane
    isOnPlane = (bsxfun(@minus, x0(i,:), x0(i+1:end,:))*n(i,:)') < 50*eps;
    if any(isOnPlane)
        % test normals are parallel
        isOnPlane = [false(i,1);isOnPlane];
        r(i) = any(abs(n(isOnPlane,:)*(n(i,:)'))> 1 - 50*eps);
    end
end
end


function newPts = mirror(pts, x0, n)
% Mirrors pts on plane
d = -2*bsxfun(@minus, pts, x0)*n';
if any(d<eps)
  warning('Found site almost on, or outside boundary! Can not guarantee for results')
end
newPts = pts + d*n;
end
