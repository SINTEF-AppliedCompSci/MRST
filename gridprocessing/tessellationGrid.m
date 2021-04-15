function G = tessellationGrid(p, t)
%Construct valid grid definition from points and tessellation list
%
% SYNOPSIS:
%   G = tessellationGrid(P, T)
%
% PARAMETERS:
%   P     - Node coordinates.  Must be an m-by-2 matrix, one row for each
%           node/point.
%
%   T     - Tessellation list: an n-by-k matrix where each row holds node
%           numbers for a k-polygon, or a n-by-1 cell array where each
%           entry holds an array of node numbers, where the array length
%           can vary for each entry
%
% RETURNS:
%   G     - Valid grid definition.
%
% EXAMPLE:
%
%   % Construct a standard Cartesian grid
%   [nx,ny] = deal(15,10);
%   [x,y] = meshgrid(linspace(0,1,nx+1),linspace(0,1,ny+1));
%   p = [x(:) y(:)];
%   n = (nx+1)*(ny+1);
%   I = reshape(1:n,ny+1,nx+1);
%   T = [
%       reshape(I(1:end-1,1:end-1),[],1)';
%       reshape(I(1:end-1,2:end  ),[],1)';
%       reshape(I(2:end,  2:end  ),[],1)';
%       reshape(I(2:end,  1:end-1),[],1)'
%       ]';
%   G = tessellationGrid(p, T);
%   clf, plotGrid(G);
%
%   % Construct a symmetric pattern of 6-point polygons
%   [dx, dy, dPhi] = deal(cos(pi/3),sin(pi/3), pi*15/180);
%   v = pi/180*[0 120 240]';
%   dv = [cos(v-dPhi) sin(v-dPhi) cos(v+dPhi) sin(v+dPhi)]/2;
%   P = [ 0           0           0           0;
%         0+dv(1,1)   0+dv(1,2)   0+dv(1,3)   0+dv(1,4);
%         1           0           1           0;
%         1+dv(2,1)   0+dv(2,2)   1+dv(2,3)   0+dv(2,4);
%         dx          dy          dx          dy;
%         dx+dv(3,1)  dy+dv(3,2)  dx+dv(3,3)  dy+dv(3,4)];
%   P1 = P(:,1:2);
%   P2 = [P([1 6:-1:2],3) -P([1 6:-1:2],4)];
%   T  = reshape(1:24,6,4)';
%
%   [p,t,n] = deal([],[],0);
%   for j=0:2
%       for i=0:4
%           p = [p; bsxfun(@plus,P1,[i 2*j*dy])];
%           p = [p; bsxfun(@plus,P2,[i-dx (2*j+1)*dy])];
%           p = [p; bsxfun(@plus,P1,[i-dx (2*j+1)*dy])];
%           p = [p; bsxfun(@plus,P2,[i 2*(j+1)*dy])];
%           t = [t; T+n]; n=n+24;
%      end
%   end
%
%   [p,ia,ic] = unique(round(p*1e5)/1e5,'rows');
%   G = tessellationGrid(p, ic(t));
%   i=repmat((1:2)',G.cells.num/2,1);
%   plotCellData(G,i(:));
%   plotFaces(G,find(any(G.faces.neighbors==0,2)),'EdgeColor','r','LineWidth',2);
%   axis tight off;
%
%   % Construct a polyhedral grid
%   p = [0 0; 1 0 ; 0.4 0.6; 1.5 0.5; 0 1; 1 1];
%   T = {[1 3 5],[3 6 5],[1 2 4 6 3]}';
%   G = tessellationGrid(p, T);
%   clf, plotGrid(G);
%
% SEE ALSO:
%   `triangleGrid`, `tetrahedralGrid`, `grid_structure`

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


   assert (size(p, 2) == 2, ...
          ['Function ''%s'' is only supported in two ', ...
           'space dimensions.'], mfilename);

   assert (nargin == 2,...
       'Must supply both a point set and a tessellation');

   if iscell(t)
      if size(t, 1) < size(t, 2)
         assert(size(t, 1) == 1, ...
            'Tesselation cell array must be one dimensional')
         t = t';
      end
      n = cellfun(@numel, t);
      tmp = arrayfun(@(n) rldecode([1:n 1]',[1 2*ones(1,n-1) 1]'), n, ...
         'uniformOutput', false);
      tt = cellfun(@(x, ix) reshape(x(ix),2,[]), ...
         t, tmp, 'UniformOutput', false);
      tt = [tt{:}];
      nc = numel(t);
   else
      n   = size(t,2);
      idx = rldecode([1:n 1]',[1 2*ones(1,n-1) 1]');
      tt  = reshape(t(:, idx)', 2, []);
      nc = size(t,1);
      n = repmat(n, [nc, 1]);
   end
   
   assert (all(max(tt) <= size(p, 1)), ...
       'Tessellation list ''T'' references invalid points.');

   assert (all(min(tt) > 0), ...
       'Tessellation list ''T'' references invalid points.');

   [fn, i]           = sort(tt);
   [fn, cf, cf]      = unique(fn', 'rows');  %#ok

   G.faces.nodes     = reshape(fn', [], 1);
   G.cells.faces     = cf;

   G.nodes.coords    = p;
   G.nodes.num       = size(p, 1);

   G.cells.num       = nc;
   G.cells.facePos   = cumsum([1; n]);

   cellNo            = rldecode(1:G.cells.num, diff(G.cells.facePos), 2).';

   G.faces.num       = double(max(G.cells.faces));
   G.faces.neighbors = accumarray([G.cells.faces, i(1,:)'], cellNo, ...
                                  [G.faces.num, 2]);

   G.faces.nodePos   = cumsum([1; repmat(2, [G.faces.num, 1])]);

   G.type            = { mfilename };
   G.griddim         = 2;

   % Uniquify nodes
   h = zeros([G.nodes.num, 1]);
   h(G.faces.nodes) = 1;

   if sum(h) ~= G.nodes.num,
      [G.nodes.coords, a, map] = unique(G.nodes.coords, 'rows');  %#ok
      G.nodes.num = size(G.nodes.coords, 1);
      G.faces.nodes = map(G.faces.nodes);
   end
end
