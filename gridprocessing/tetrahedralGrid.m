function G = tetrahedralGrid(p, varargin)
%Construct valid grid definition from points and tetrahedron list
%
% SYNOPSIS:
%   G = tetrahedralGrid(P)
%   G = tetrahedralGrid(P, T)
%
% PARAMETERS:
%   P - Node coordinates.  Must be an m-by-3 matrix, one row for each
%       node/point.
%
%   T - Tetrahedron list (point tesselation): an n-by-4 matrix where each
%       row holds node numbers for a tetrahedron.
%
%       OPTIONAL.  Default value: T = DELAUNAY3(P(:,1), P(:,2), P(:,3))
%
% RETURNS:
%   G  - Valid grid definition.
%
% SEE ALSO:
%   `triangleGrid`, `grid_structure`

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


   assert (size(p, 2) == 3, ...
          ['Function ''%s'' is only supported in three ', ...
           'space dimensions.'], mfilename);

   t = tesselate(p, varargin{:});

   % This must be implemented to make code robust wrt t.
   t = permuteNodes(t, p);

   [fn,i]            = sort(reshape(t(:,[2,1,3, 1,2,4, 3,1,4, 2,3,4])', 3, []));
   [fn, cf, cf]      = unique(fn', 'rows');
   G.faces.nodes     = reshape(fn', [], 1);
   G.cells.faces     = cf;

   G.nodes.coords    = p;
   G.nodes.num       = size(p, 1);

   G.cells.num       = size(t, 1);
   G.cells.indexMap  = (1:G.cells.num) .';
   G.cells.facePos   = cumsum([1;repmat(4, [G.cells.num, 1])]);
   cellNo            = rldecode(1:G.cells.num, diff(G.cells.facePos), 2).';


   G.faces.num       = double(max(G.cells.faces));
   G.faces.neighbors = accumarray([G.cells.faces,  ...
                       1+(any(diff(i)'==1, 2))], cellNo, ...
                       [G.faces.num, 2]);
   G.faces.nodePos   = cumsum([1;repmat(3, [G.faces.num, 1])]);
   G.type            = { mfilename };
   G.griddim         = 3;

   % Uniquify nodes
   h=zeros(G.nodes.num, 1);
   h(G.faces.nodes) = 1;
   if sum(h) ~= G.nodes.num,
      [G.nodes.coords, a, map] = unique(G.nodes.coords, 'rows');
      G.nodes.num = size(G.nodes.coords, 1);
      G.faces.nodes = map(G.faces.nodes);
   end
end

%--------------------------------------------------------------------------

function t = tesselate(p, varargin)
   if nargin < 2,
      if exist('delaunay3', 'file')
         t = delaunay3(p(:,1), p(:,2), p(:,3)); %#ok, we check for this
      else
         t = delaunay(p(:,1), p(:,2), p(:,3));
      end

   elseif isnumeric(varargin{1}) && (size(varargin{1}, 2) ==4),
      t = varargin{1};

      assert (all(max(t) <= size(p, 1)), ...
              'Tetrahedron list ''T'' references invalid points.');

      assert (all(min(t) > 0), ...
              'Tetrahedron list ''T'' references invalid points.');

   else
      error(['Input parameter ''T'' does not appear to be a valid ', ...
             'tesselation of the points ''P''.']);
   end
end

%--------------------------------------------------------------------------

function t = permuteNodes(t, p)
%T = permuteNodes(T,P)
%
%  permute nodes in tetrahedra such that the triple product
%  p(t(:,2:4),:)-p(t(:,1),:) (i.e., the volume) is positive.

   v = tripleProduct(t, p);
   t(v>0,1:2) = t(v>0, [2,1]);
end

%--------------------------------------------------------------------------

function v = tripleProduct(t, p)
%
% Compute triple product (v2 x v3)·v4 where vi = p(t(:,i),:)-p(t(:,1),:).

   x = p(:,1); y = p(:,2); z = p(:,3);
   X = x(t);   Y = y(t);   Z = z(t);

   dx = bsxfun(@minus, X(:, 2:4), X(:,1));
   dy = bsxfun(@minus, Y(:, 2:4), Y(:,1));
   dz = bsxfun(@minus, Z(:, 2:4), Z(:,1));

   v= dot(cross([dx(:,1), dy(:,1), dz(:,1)], ...
             [dx(:,2), dy(:,2), dz(:,2)]), ...
             [dx(:,3), dy(:,3), dz(:,3)], 2);
end
