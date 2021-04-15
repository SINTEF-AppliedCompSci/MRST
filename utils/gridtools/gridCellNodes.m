function [n, pos] = gridCellNodes(G, c, varargin)
%Extract nodes per cell in a particular set of cells
%
% SYNOPSIS:
%   [n, pos] = gridCellNodes(G, c)
%   [n, pos] = gridCellNodes(G, c, 'pn1', pv1, ...)
%
% PARAMETERS:
%   G       - Grid structure.
%
%   c       - Cells for which to extract unique nodes.  Array of numeric
%             cell indices.
%
% OPTIONAL PARAMETERS:
%   'unique' - Whether or not to compute *unique* nodes per cell. Boolean
%             flag (Logical).  Default value: unique = true (do compute
%             uniuqe nodes without repetitions per cell).
%
% NOTE:
%   Computing unique nodes (i.e., setting option 'unique' to `true`)
%   increases the computational complexity of function `gridCellNodes`
%   because this leads to invoking function `sortrows`.
%
% RETURNS:
%   n   - Nodes per cell.  Unique nodes if option `unique` is `true`.
%   pos - Indirection map into n.  The nodes of cell `c(i)` are in positions
%         `p(i) : p(i+1) - 1`.
%
% EXAMPLE:
%   g = cartGrid([2, 2, 2])
%   [u, pu] = gridCellNodes(g, 1)
%   [n, pn] = gridCellNodes(g, 1, 'unique', false)
%
% SEE ALSO:
%   `sortrows`.

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


    opt = struct('unique', true);
    opt = merge_options(opt, varargin{:});

    % Find number of faces per cell
    nf = diff([G.cells.facePos(c), G.cells.facePos(c+1)], [],2);

    % Find the cell index of each face
    cellno = rldecode(1:numel(c), nf, 2) .';

    % Number of nodes per face
    nnode = diff(G.faces.nodePos);

    % Find faces of cell subset
    cf = G.cells.faces(mcolon(G.cells.facePos( c ), ...
                              G.cells.facePos(c+1) - 1),1);
    % Find node indices of the faces of the cell subset through indrection
    % map.
    ni = mcolon(G.faces.nodePos(cf), ...
                G.faces.nodePos(cf+1)-1)';

    W = [rldecode(cellno, nnode(cf)), G.faces.nodes(ni)];
    if opt.unique,
       % Caller requested a set of unique nodes per cell.
       % Incur (a usually) expensive SORTROWS call.
       W = rlencode(sortrows(W));
    end

    pos = cumsum([1; accumarray(W(:,1),1)]);
    n = W(:,2);

    assert (numel(pos) == numel(c) + 1);
end
