function Gf = gridFracturePlane_n(fracplane, varargin)
% Grids a fracture plane defined by a set of coplanar points. The coplanar
% points are assumed to define a convex polygon.

opt = struct('type', 'triangle', 'useDistmesh', false, 'elementSize', 1);
opt = merge_options(opt, varargin{:});
%
flag = 1*strcmp(opt.type,'triangle') + 2*strcmp(opt.type,'pebi');
%
xyp = rotatePlane(fracplane.points,[0 0 1]);
z = xyp(1,3);

xyp = xyp(:,1:2);
% xyp = gridPolygon2D(xyp);
%
if opt.useDistmesh
    endp_xy = xyp(1:size(fracplane.points,1),1:2);
%     h = min(sqrt(sum(diff([endp_xy;endp_xy(1,:)],1).^2,2))); % min side length
    h = sqrt(opt.elementSize);
    expand = h;
    esize = h/2;
    bbox = [min(endp_xy(:,1)) - expand, min(endp_xy(:,2)) - expand; ...
        max(endp_xy(:,1)) + expand, max(endp_xy(:,2)) + expand];
    try
        [p,t] = distmesh_2d(@dpoly, @huniform, esize, bbox, 500, ...
            endp_xy(:,1:2),endp_xy(:,1:2));
%         close gcf
        Gf = triangleGrid(p,t);
        Gf = computeGeometry(Gf);
    catch
        Gf = triangleGrid(xyp);
    end
else
    Gf = triangleGrid(xyp);
    Gf = computeGeometry(Gf);
end

if flag == 2
    try
        Gf = pebi(Gf);
    catch
        fprintf('\nPointset not favourable for PEBI Grids. Proceeding with triangle grid.\n');
    end
end

nc = [Gf.nodes.coords,repmat(z,Gf.nodes.num,1)];
Gf = makeLayeredGrid(Gf,1);
Gf = computeGeometry(Gf);
%
unitNormal = fracplane.normal/norm(fracplane.normal);
tdist = fracplane.aperture/2;

use_points = rotatePlane(nc,+unitNormal);
if any(use_points(:,3)<0)
    use_points = rotatePlane(nc,-unitNormal);
end
% Projecting on plane, some numerical error will always be there
% pdash = [fracplane.points;fracplane.points(1,:)];
% basisCoefficients = bsxfun(@minus,use_points,pdash(1,:))*null(unitNormal);
% use_points = bsxfun(@plus,basisCoefficients*null(unitNormal).', pdash(1,:));
%
% use_points(:,3) = use_points(:,3) + abs(min(use_points(:,3)));
pminus = use_points - repmat(unitNormal,size(use_points,1),1)*tdist;
pplus = use_points + repmat(unitNormal,size(use_points,1),1)*tdist;


Gf.nodes.coords = [pminus;pplus];
Gf = computeGeometry(Gf);

return