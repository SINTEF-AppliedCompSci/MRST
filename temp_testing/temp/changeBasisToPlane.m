function points_xy = changeBasisToPlane(points,varargin)

opt = struct('origin', [], 'basis', []);
opt = merge_options(opt, varargin{:});
if isempty(opt.origin)
    opt.origin = points(1,:);
else
    [~,loc] = ismember(opt.origin,points,'rows');
    if loc~=1 && loc~=0
        points = [points(loc:end,:);points(1:loc-1,:)];
    end
end
newPosVec = points - repmat(opt.origin,size(points,1),1);
% Gram Schmidt orthonormalisation
if isempty(opt.basis)
    basis1 = (newPosVec(2,:))/norm(newPosVec(2,:));
    basis2 = newPosVec(3,:) - basis1*dot(newPosVec(3,:),basis1);
    basis2 = basis2/norm(basis2);
    basis3 = cross(basis1,basis2);
else
    basis1 = opt.basis(1,:);
    basis2 = opt.basis(2,:);
    basis3 = opt.basis(3,:);
end
%
points_xy = zeros(size(points));
for i = 1:size(points,1)
    points_xy(i,:) = [dot(newPosVec(i,:),basis1),...
                      dot(newPosVec(i,:),basis2),...
                      dot(newPosVec(i,:),basis3)];
end

d1 = -points(1,:);
d1d = [dot(d1,basis1),dot(d1,basis2),dot(d1,basis3)];
assert(abs(norm(d1) - norm(d1d))<eps*100,'Check coplanarity of points.'); % norms must be equal

return