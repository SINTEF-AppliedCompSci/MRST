function dz = getDZ(G, c)
% Compute DZ of cell c in grid G. The face direction indicator in
% cell.faces(:,2) is required
%
% EXAMPLE:
%   G = cartGrid([10, 10, 5], [100, 100, 10]);
%   G = computeGeometry(G);
%   dz = getDZ(G, 1);

    fPos    = G.cells.facePos(c) : G.cells.facePos(c+1)-1;
    faces   = G.cells.faces(fPos,1);
    faceDir = G.cells.faces(fPos,2);
    fCenter6 = G.faces.centroids(faces(faceDir == 6), :);
    fCenter5 = G.faces.centroids(faces(faceDir == 5), :);
    dz = norm(fCenter6 - fCenter5, 2);
end