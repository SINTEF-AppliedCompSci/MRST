function fracplanes = getPlaneNormals(fracplanes,varargin)

opt = struct('normalize',false);
opt = merge_options(opt, varargin{:});

for i = 1:numel(fracplanes)
    points = [fracplanes(i).points;fracplanes(i).points(1,:)];
    
    diffp = diff(points,1);
    normal = cross(diffp(1,:), diffp(2,:));
    if opt.normalize
        normal = normal/norm(normal);
    end
    
    nBoundPlanes = size(points,1)-1;
    boundNormals = zeros(nBoundPlanes,3);
    for j = 1:nBoundPlanes
        boundNormals(j,:) = cross(-diffp(j,:),normal);
        if opt.normalize
            boundNormals(j,:) = boundNormals(j,:)/norm(boundNormals(j,:));
        end
    end
    fracplanes(i).normal = normal;
    fracplanes(i).boundNormals = boundNormals;
end