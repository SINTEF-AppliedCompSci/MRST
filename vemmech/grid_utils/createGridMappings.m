function w  = createGridMappings(G)
%
%
% SYNOPSIS:
%   function w  = createGridMappings(g)
%
% DESCRIPTION: Create preliminary structures which is used to conveniently
% set up the grid mapping, see full_grid_structure.
%
% PARAMETERS:
%   G - Grid structure
%
% RETURNS:
%   w - preliminary structure to be used in createAugmentedGrid
%
% EXAMPLE:
%
% SEE ALSO:
%
    
    cellno   = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
    col      = 1 + (cellno == G.faces.neighbors(G.cells.faces(:,1), 2));
    nhfaces  = G.cells.facePos(end)-1;
    hfaces   = accumarray([G.cells.faces(:,1), col], 1:nhfaces);
    hfaces   = rldecode(hfaces, diff(G.faces.nodePos));
    cells    = rldecode(G.faces.neighbors, diff(G.faces.nodePos));
    nodes    = repmat(G.faces.nodes, [2,1]);
    faces    = repmat(rldecode(1:G.faces.num, diff(G.faces.nodePos),2)', [2,1]);
    i        = cells~=0;
    w        = [cells(i), nodes(i), hfaces(i), faces(i)];
    w        = double(sortrows(w));
    
end
