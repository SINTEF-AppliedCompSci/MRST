function [P, faceNo, cellNo, w] = getFaceFromCellInterpolator(G)
    % Get nodes
    faceNo = G.cells.faces(:, 1);
    pos = G.cells.facePos;
    cellNo = rldecode((1:G.cells.num)', diff(pos));
    
    % One over distance weighting
    w = 1./sqrt(sum((G.cells.centroids(cellNo, :) - G.faces.centroids(faceNo, :)).^2, 2));
    
    % Divide by sum of weights for partition of unity
    sumw = accumarray(faceNo, w);
    w = w./sumw(faceNo);
    
    % n_n by n_c matrix constructing values on nodes by cells
    P = sparse(faceNo, cellNo, w, G.faces.num, G.cells.num);
end