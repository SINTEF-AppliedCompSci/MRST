function [P, nodeNo, cellNo, w] = getNodeFromCellInterpolator(G)
    % Get nodes
    [nodeNo, pos] = gridCellNodes(G, 1:G.cells.num, 'unique', true);
    cellNo = rldecode((1:G.cells.num)', diff(pos));
    
    % One over distance weighting
    w = 1./sqrt(sum((G.cells.centroids(cellNo, :) - G.nodes.coords(nodeNo, :)).^2, 2));
    
    % Divide by sum of weights for partition of unity
    sumw = accumarray(nodeNo, w);
    w = w./sumw(nodeNo);
    
    % n_n by n_c matrix constructing values on nodes by cells
    P = sparse(nodeNo, cellNo, w, G.nodes.num, G.cells.num);
end