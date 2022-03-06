function [vec, p, f] = getDirectionToClosestFace(G, pnt, faceIx)
% Undocumented utility function for trajectory optimization

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
c = G.faces.centroids(faceIx,:);
d  = bsxfun(@minus, c, pnt);
nds = dot(d,d,2);
[~, ix] = min(nds);
vec = d(ix,:);
f = faceIx(ix);
p = c(ix,:);
end