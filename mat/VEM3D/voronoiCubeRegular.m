function G = voronoiCubeRegular(gridDim,gridLim,perturb)

    boundary = [0,          0         , 0         ; ...
                gridLim(1), 0         , 0         ; ...
                gridLim(1), gridLim(2), 0         ; ...
                0         , gridLim(2), 0         ; ...
                0         , 0         , gridLim(3); ...
                gridLim(1), 0         , gridLim(3); ...
                gridLim(1), gridLim(2), gridLim(3); ...
                0         , gridLim(2), gridLim(3)];
            
    dx = gridLim(1)/gridDim(1);
    dy = gridLim(2)/gridDim(2);
    dz = gridLim(3)/gridDim(3);
            
    G = cartGrid(gridDim-1, gridLim-[dx,dy,dz]);
    

    
    X = G.nodes.coords;
%     X = X(X(:,1) ~= 0 & X(:,1) ~= gridLim(1) & ...
%           X(:,2) ~= 0 & X(:,2) ~= gridLim(2) & ...
%           X(:,3) ~= 0 & X(:,3) ~= gridLim(3),:);
%     
    
    X = bsxfun(@plus, X, [dx, dy, dz]/2);
    n = size(X,1);
% 
%     xEnd = X(:,1) == 0 | X(:,1) == gridLim(1);
%     nX = sum(~xEnd);
%     X(~xEnd,1) = X(~xEnd,1) + (rand(nX,1)*dx-dx/2)*perturb;
    X(:,1) = X(:,1) + (rand(n,1)*dx-dx/2)*perturb;
    
    
%     yEnd = X(:,2) == 0 | X(:,2) == gridLim(2);
%     nY = sum(~yEnd);
%     X(~yEnd,2) = X(~yEnd,2) + (rand(nX,1)*dy-dy/2)*perturb;
    X(:,2) = X(:,2) + (rand(n,1)*dy-dy/2)*perturb;
    
%     zEnd = X(:,3) == 0 | X(:,3) == gridLim(3);
%     nZ = sum(~zEnd);
%     X(~zEnd,3) = X(~zEnd,3) + (rand(nZ,1)*dz-dz/2)*perturb;
    X(:,3) = X(:,3) + (rand(n,1)*dz-dz/2)*perturb;
    
    G = voronoi3D(X, boundary);

end

