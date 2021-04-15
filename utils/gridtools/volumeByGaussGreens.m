function V = volumeByGaussGreens(G)
%Compute cell volume by means of Gauss-Greens' formula
%
% SYNOPSIS:
%   V = volumeByGaussGreens(G)
%
% PARAMETERS:
%   G - Grid structure.
%
% RETURNS:
%   V - G.cells.num-by-1 array of cell volumes computed by means of
%       Gauss-Greens formula.  This is a useful consistency check when
%       developing new grid processing tools.  In most cases, this value
%       should not deviate too much from the 'G.cells.volumes' value
%       computed by the 'computeGeometry' function.
%
% SEE ALSO:
%   `computeGeometry`.

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


   cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';

   sgn = 1-2*(G.faces.neighbors(G.cells.faces(:,1),2)==cellNo);
   R   = G.faces.centroids(G.cells.faces(:,1),:);% - G.cells.centroids(cellNo,:);
   N   = bsxfun(@times, G.faces.normals(G.cells.faces(:,1),:), sgn);

   %% Reshape N and R to a block-diagonal matrices.
   d   = G.griddim;
   [i,j] = blockDiagIndex(repmat(d, [G.cells.num, 1]), diff(G.cells.facePos));

   % reshape
   R  = sparse(i,j,R')';
   N  = sparse(i,j,N')';

   %% Gauss-Greens formula
   % Apply Gauss-Greens formula to vector field [x,0,0],
   %
   %     v = \int_C1 div([x,0,0]) dv = \int_{dC1} x·n ds.
   %
   % For polygons and polyhedra we have:
   %
   %     \int_{dC1} x·n ds = R(:,1)'*N(:,1).
   %
   % The same holds for [0,y,0] and [0,0,z].

   V = reshape(spdiags(R'*N, 0), d, [])';
end
