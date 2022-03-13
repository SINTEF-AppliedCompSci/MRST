function h = plotCellNumbers(g, varargin)
% Debug function which plots cell numbers on a (subset) of cells

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

   if mod(nargin, 2),
      c = 1:g.cells.num;
   else
      c = varargin{1};
      varargin = varargin(2:end);
   end
   if isfield(g, 'parent'),
      dispif(mrstVerbose, 'Adding (some) geometry to coarse grid.');
      if ~isfield(g.cells, 'centroids'),
            g = computeCoarseGeometry(g);
      end
   elseif ~isfield(g.cells, 'centroids'),
      dispif(mrstVerbose, 'Adding geometry to grid.');
      g = computeGeometry(g);
   else
      error('');
   end
   v = ishold;hold on;
   h=text(g.cells.centroids(c,1), g.cells.centroids(c,2), num2str(c(:)), varargin{:});
   if ~v, hold off; end
end


function cg = computeCoarseGeometry(cg)
   if ~isfield(cg.parent.faces, 'centroids'),
      dispif(mrstVerbose, 'Adding geometry to fine grid.');
      cg.parent = computeGeometry(cg.parent);
   end

   cg.faces.centroids = zeros(cg.faces.num, cg.griddim);
   faceno = rldecode(1:cg.faces.num, diff(cg.faces.connPos), 2)';
   for i=1:cg.griddim,
      cg.faces.centriods(:,i) = accumarray(faceno, cg.parent.faces.centroids(cg.faces.fconn, i));
   end
   cg.faces.centriods = bsxfun(@rdivide, cg.faces.centriods, accumarray(faceno, cg.parent.faces.areas(cg.faces.fconn)));

   cg.cells.centroids = zeros(cg.cells.num, cg.griddim);
   for i=1:cg.griddim,
      cg.cells.centroids(:,i) = accumarray(cg.partition, cg.parent.cells.centroids(:, i).*cg.parent.cells.volumes);
   end
   cg.cells.centroids = bsxfun(@rdivide, cg.cells.centroids, accumarray(cg.partition, cg.parent.cells.volumes));
end
