function Gf = gridFracturePlane(fracplane, nodes, varargin)
% Grids a fracture plane defined by a set of coplanar points. The coplanar
% points are assumed to define a convex polygon.

opt = struct('type', 'triangle');
opt = merge_options(opt, varargin{:});

flag = 1*strcmp(opt.type,'triangle') + 2*strcmp(opt.type,'pebi');
%
endp = fracplane.points;
endp_xy = rotatePlane(endp,[0 0 1]); % End points of rotated plane // to x-y plane
z = endp_xy(1,3);
assert(numel(uniquetol(endp_xy(:,3),eps*100,'DataScale', 1))==1,...
    'Incorrect axis rotation');
%
use_points = unique([endp;nodes],'rows','stable');
points_xy = rotatePlane(use_points,[0 0 1]); % Make plane // to x-y plane
%

% points_xy = round([endp_xy;points_xy],2);
points_xy = unique(points_xy(:,1:2),'rows','stable');
% 
endp_xy = endp_xy(:,1:2);
h = min(sqrt(sum(diff([endp_xy;endp_xy(1,:)],1).^2,2))); % min side length
expand = h/4;
esize = h/4;
bbox = [min(endp_xy(:,1)) - expand, min(endp_xy(:,2)) - expand; ...
        max(endp_xy(:,1)) + expand, max(endp_xy(:,2)) + expand];
[p,t] = distmesh_2d(@dpoly, @huniform, esize, bbox, 1000, ...
        endp_xy(:,1:2),endp_xy(:,1:2));
[p,t] = fixmesh(p,t);
close gcf
%
Gf = triangleGrid(p,t);
if flag == 2
    Gf = pebi(Gf);
end
nc = [Gf.nodes.coords,repmat(z,Gf.nodes.num,1)];
Gf = makeLayeredGrid(Gf,1);
%
unitNormal = fracplane.normal/norm(fracplane.normal);
tdist = fracplane.aperture/2;

use_points = rotatePlane(nc,+unitNormal);
use_points(:,3) = use_points(:,3) + abs(min(use_points(:,3)));


pminus = use_points - repmat(unitNormal,size(use_points,1),1)*tdist;
pplus = use_points + repmat(unitNormal,size(use_points,1),1)*tdist;


Gf.nodes.coords = [pminus;pplus];
Gf = computeGeometry(Gf);

return