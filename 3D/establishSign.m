function Sign = establishSign(G, fracplanes, varargin)

opt = struct('tolerance', eps*100);
opt = merge_options(opt, varargin{:});
tol = opt.tolerance;

Sign = struct;
fcents = G.faces.centroids;
ncoords = G.nodes.coords;
for i = 1:numel(fracplanes)
    normal = fracplanes(i).normal;
    points = [fracplanes(i).points;fracplanes(i).points(1,:)];
    p = mean(points,1);
    A = normal(1); B = normal(2); C = normal(3);
    D = -dot(normal,p);
    
    sign = A*fcents(:,1) + B*fcents(:,2) + C*fcents(:,3) + D;
    sign(abs(sign)>tol) = sign(abs(sign)>tol)./abs(sign(abs(sign)>tol));
    sign(isnan(sign)) = 0;
    sign(abs(sign)<tol) = 0;
    Sign(i).FaceSign(:,1) = sign;
    
    sign = A*ncoords(:,1) + B*ncoords(:,2) + C*ncoords(:,3) + D;
    sign(abs(sign)>tol) = sign(abs(sign)>tol)./abs(sign(abs(sign)>tol));
    sign(isnan(sign)) = 0;
    sign(abs(sign)<tol) = 0;
    Sign(i).NodeSign(:,1) = sign;
    
    for j = 1:size(fracplanes(i).boundNormals,1)
        normal = fracplanes(i).boundNormals(j,:);
        A = normal(1); B = normal(2); C = normal(3);
        D = -dot(normal,mean(points(j:j+1,:),1));
        
        sign = A*fcents(:,1) + B*fcents(:,2) + C*fcents(:,3) + D;
        sign(abs(sign)>tol) = sign(abs(sign)>tol)./abs(sign(abs(sign)>tol));
        sign(isnan(sign)) = 0;
        sign(abs(sign)<tol) = 0;
        Sign(i).FaceSign(:,j+1) = sign;
        
        sign = A*ncoords(:,1) + B*ncoords(:,2) + C*ncoords(:,3) + D;
        sign(abs(sign)>tol) = sign(abs(sign)>tol)./abs(sign(abs(sign)>tol));
        sign(isnan(sign)) = 0;
        sign(abs(sign)<tol) = 0;
        Sign(i).NodeSign(:,j+1) = sign;
    end
end

return