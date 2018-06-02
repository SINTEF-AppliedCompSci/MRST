function [dx, dy, dz] = cellDims(G, ix)
% cellDims -- Compute physical dimensions of all cells in single well
%
% SYNOPSIS:
%   [dx, dy, dz] = cellDims(G, ix)
%
% PARAMETERS:
%   G  - Grid data structure.
%   ix - Cells for which to compute the physical dimensions
%
% RETURNS:
%   dx, dy, dz -- [dx(k) dy(k)] is bounding box in xy-plane, while dz(k) =
%                 V(k)/dx(k)*dy(k)

n = numel(ix);
[dx, dy, dz] = deal(zeros([n, 1]));

ixc = G.cells.facePos;
ixf = G.faces.nodePos;

for k = 1 : n
   c = ix(k);                                     % Current cell
   f = G.cells.faces(ixc(c) : ixc(c + 1) - 1, 1); % Faces on cell
   e = mcolon(ixf(f), ixf(f + 1) - 1);            % Edges on cell

   nodes  = unique(G.faces.nodes(e, 1));          % Unique nodes...
   coords = G.nodes.coords(nodes,:);            % ... and coordinates

   % Compute bounding box
   m = min(coords);
   M = max(coords);

   % Size of bounding box
   dx(k) = M(1) - m(1);
   if size(G.nodes.coords, 2) > 1
      dy(k) = M(2) - m(2);
   else
      dy(k) = 1;
   end

   if size(G.nodes.coords, 2) > 2
      dz(k) = G.cells.volumes(ix(k))/(dx(k)*dy(k));
   else
      dz(k) = 1;
   end
end
