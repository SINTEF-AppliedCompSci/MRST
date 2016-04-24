function p = rotatePlane(points, normal, varargin)
% Rotates plane defined by points to align with plane defined by input
% argument normal. Rotation is along the axis defined by the intersection
% of the two planes. Returns points as it is if
% both planes are parallel.

opt = struct('useNormal',[]);
opt = merge_options(opt,varargin{:});
assert(size(points,1)>=3,'3 or more points needed to define a plane!');
% assert(all(iscoplanar(points)),'All Points passed must be coplanar');
if isempty(opt.useNormal)
    diffp = diff(points,1);
    normal2 = cross(diffp(1,:), diffp(2,:));
    normal2 = normal2/norm(normal2); 
else
    normal2 = opt.useNormal/norm(opt.useNormal);
end
normal = normal/norm(normal);
if isequal(abs(normal),abs(normal2))
    p = points;
    return
end

rotAngle = -acos(dot(normal2,normal)); %rotation angle
rotAxis = -cross(normal,normal2); %rotation axis

M = makehgtform('axisrotate', rotAxis, rotAngle);
p = points*M(1:3,1:3);

return