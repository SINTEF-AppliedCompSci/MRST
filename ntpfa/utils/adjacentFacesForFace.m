function faces = adjacentFacesForFace(G, f, cts)
%Undocumented Utility Function

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

    if nargin < 3
        cts = 2;
    end
    active = false(G.nodes.num, 1);
    active(gridFaceNodes(G, f)) = true;

    % Count number of nodes shared per face
    nPos = rldecode(1:G.faces.num, diff(G.faces.nodePos), 2) .';
    counts = accumarray(nPos, active(G.faces.nodes(:, 1)));
    % Take as active anyone with more than cts overlapping nodes
    faces = find(counts >= cts);
    faces = faces(faces~=f);
end
