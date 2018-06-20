function [dx, dy, dz] = cellDimsCG(G,ix)
% cellDims -- Compute physical dimensions of all cells in single well
%
% SYNOPSIS:
%   [dx, dy, dz] = cellDims(G, ix)
%
% PARAMETERS:
%   G  - Grid data structure.
%   ix - Cells for which to compute the physical dimensions (bounding
%        boxes).
%
% RETURNS:
%   dx, dy, dz -- Size of bounding box for each cell.  In particular,
%                 [dx(k),dy(k),dz(k)] is Cartesian BB for cell ix(k).
n = numel(ix);
[dx, dy, dz] = deal(zeros([n, 1]));
ixc = G.cells.facePos;
dim = size(G.faces.centroids, 2);

for k = 1 : n
    c = ix(k);                                     % Current cell
    f = G.cells.faces(ixc(c) : ixc(c + 1) - 1, 1); % Faces on cell
    fc = G.faces.centroids(f, :);
    delta = max(fc) - min(fc);
    dx(k) = delta(1);
    dy(k) = delta(2);
    if dim > 2
        dz(k) = delta(3);
    else
        dz(k) = 1;
    end
end
