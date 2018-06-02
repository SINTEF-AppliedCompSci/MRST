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

for k = 1 : n
    c = ix(k);                                     % Current cell
    f = G.cells.faces(ixc(c) : ixc(c + 1) - 1, 1); % Faces on cell
    assert(numel(f)==6);
    [~, ff] = sortrows(abs(G.faces.normals(f,:)));
    f = f(ff(end:-1:1));
    dx(k) = 2*G.cells.volumes(c) / (sum(G.faces.areas(f(1:2))));
    dy(k) = 2*G.cells.volumes(c) / (sum(G.faces.areas(f(3:4))));
    dz(k) = 2*G.cells.volumes(c) / (sum(G.faces.areas(f(5:6))));
end
