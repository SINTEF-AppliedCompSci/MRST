function [n, pos] = gridFaceNodes(G, f)
% Find nodes corresponding to a set of faces
%
% SYNOPSIS:
%   [n, pos] = gridFaceNodes(G, c)
%
% PARAMETERS:
%   G    - Grid structure
%   c    - Cells where the fine nodes are desired
%
%
% RETURNS:
%   pos  - indirection map into `n`. The nodes of face `f(i)` is found at
%          positions `p(i):p(i+1)-1` in `n`
%   n    - node positions in `G`

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


    % Number of nodes per face
    nnode = diff(G.faces.nodePos);

    ni = mcolon(G.faces.nodePos(f), ...
                G.faces.nodePos(f+1)-1)';

    pos = cumsum([1; double(nnode(f))]);
    n = G.faces.nodes(ni);
end
