function w  = createGridMappings(G)
% Add preliminary mappings to be used in `createAugmentedGrid`
%
% SYNOPSIS:
%   function w  = createGridMappings(g)
%
% DESCRIPTION: 
% Create preliminary structures which is used to conveniently set up the
% grid mapping, see `extended_grid_structure`.
%
% PARAMETERS:
%   G - Grid structure
%
% RETURNS:
%   w - preliminary structure to be used in `createAugmentedGrid`
%
% SEE ALSO:
%   `createAugmentedGrid`.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

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
