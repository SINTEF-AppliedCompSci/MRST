function [isit, disit] = iscoplanar(pxyz,varargin)

opt = struct('normal', [], 'point', [], 'distTolerance', eps);
opt = merge_options(opt, varargin{:});
dtol = opt.distTolerance;
normal = opt.normal;
p = opt.point;
[x, y, z] = deal(pxyz(:,1),pxyz(:,2),pxyz(:,3));

% Test 1: Very Strict

warning('off','MATLAB:triangulation:EmptyTri3DWarnId');
tri = delaunayTriangulation(x,y,z);
if isempty(tri.ConnectivityList);
    isit = ones(size(x));
else
    pp = [x,y,z];
    isit = ~ismember(pp,tri.Points);
end
warning('on','MATLAB:triangulation:EmptyTri3DWarnId');

% Test 2: Distance Based
if isempty(p)
    p = [x(1),y(1),z(1)];
end
if isempty(normal)
    diffp = diff([x,y,z],1);
    normal = cross(diffp(1,:), diffp(2,:));
    normal = normal/norm(normal);
end

A = normal(1); B = normal(2); C = normal(3);
D = -dot(normal,p);

dist = abs((A*x + B*y + C*z + D)./norm(normal));
disit = dist<=dtol;
isit = logical(isit);
% disit(isit) = true;
return

