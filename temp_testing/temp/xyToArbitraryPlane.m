function p = xyToArbitraryPlane(points, planep, varargin)

opt = struct('origin', planep(1,:), 'normal', [], 'basis', []);
opt = merge_options(opt, varargin{:});
normal = opt.normal;
origin = opt.origin;

assert(size(planep,1)>=3,'3 or more points needed to define a plane!');
%
if size(points,2)<3
    p = [points,zeros(size(points,1),1)];
else
    p = points;
end

newPosVec = p + repmat(origin,size(p,1),1);

if isempty(normal)
    diffp = diff(planep,1);
    normal = cross(diffp(1,:), diffp(2,:));
end

if isempty(opt.basis)
    basis1 = [1 0 0];
    basis2 = [0 1 0];
    basis3 = [0 0 1];
else
    basis1 = opt.basis(1,:);
    basis2 = opt.basis(2,:);
    basis3 = normal;
end

for i = 1:size(points,1)
    p(i,:) = [dot(newPosVec(i,:),basis1),...
        dot(newPosVec(i,:),basis2),...
        dot(newPosVec(i,:),basis3)];
end

return