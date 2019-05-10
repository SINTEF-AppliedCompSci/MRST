function N = findStreamlineNeighborshipCart(G)
% Build (n x 2*d) -array of neighbors for each cell in (Cartesian) grid G.
   cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
   col    = 1 +   (cellNo == G.faces.neighbors(G.cells.faces(:,1), 1));
   c      = G.faces.neighbors(double(G.cells.faces(:,1)) + G.faces.num* (col-1));
   N      = accumarray([cellNo, G.cells.faces(:,2)], c);
end