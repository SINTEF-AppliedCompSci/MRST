function inside = isPointInsideGrid(G, p )
%Determine whether a set of points lies within the circumference of a grid
%
% SYNOPSIS:
%    ii = isPointInsideGrid(G, p)
%
% PARAMETERS:
%    G  - grid structure
%    p  - nx2 array of points
%
% RETURNS:
%    ii - logical vector with length equal size(p,1), true if point is
%         inside and false if not
%
% Disclaimer: works only for 2D grids

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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

assert(G.griddim==2, 'Implemented only for 2D grids');
assert(size(p,2)==2, 'Point list must consist of 2D points');

% -- Extract nodes on circumference ---------------------------------------
f = boundaryFaces(G);
n = gridFaceNodes(G, f);

% -- Form a polygon from the boundary nodes -------------------------------
% The nodes come in pairs. Sort the pairs so that the second node in the
% first pair is the first node in the second pair, etc. Save only the first
% node of the pair since each node will appear twice.
m = reshape(n,2,[])';
num  = size(m,1);
i    = false(num,1); i(1) = true;
node = zeros(num,1); node(1) = m(1,1);
rows = [1 2];
for k=2:size(m,1)
   node(k) = m(i,rows(2));
   m(i,rows(2)) = 0;
   i = m(:,rows(1))==node(k);
   if ~any(i)
      rows = rows([2 1]);
      i = m(:,rows(1))==node(k);
   end
end
node(end) = node(1);

% -- Are points inside the closed polygon? --------------------------------
inside = inpolygon(p(:,1), p(:,2), ...
   G.nodes.coords(node,1),G.nodes.coords(node,2));
end
