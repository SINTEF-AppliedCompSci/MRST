function [p, n] = cellnode(G, c)
% Find nodes corresponding to a set of cells
%
% SYNOPSIS:
%   [p, n] = cellnode(G, c)
%
% PARAMETERS:
%   G    - Grid structure
%   c    - Cells where the fine nodse are desired
%
%
% RETURNS:
%   p    - indirectionmap into n. The nodes of cell c(i) is found at
%   positions p(i):p(i+1)-1 in n
%   n    - node positions in G

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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

    nf = diff([G.cells.facePos(c), G.cells.facePos(c+1)], [],2);
    cellno = rldecode(1:numel(c), nf, 2) .';
    nnode = diff(G.faces.nodePos);
    cf = G.cells.faces(mcolon(G.cells.facePos( c ), ...
                              G.cells.facePos(c+1) - 1),1);
    ni = mcolon(G.faces.nodePos(cf), ...
                G.faces.nodePos(cf+1)-1)';
    W = [rldecode(cellno, nnode(cf)), G.faces.nodes(ni)];
    W = rlencode(sortrows(W));
    p = cumsum([1; accumarray(W(:,1),1)]);
    n = W(:,2);

    assert (numel(p) == numel(c) + 1);
end
