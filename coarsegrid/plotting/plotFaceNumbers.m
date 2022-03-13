function h = plotFaceNumbers(g, varargin)
% Debug utility which plots face numbers on a (subset) of faces.

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
      f = 1:g.faces.num;
   else
      f = varargin{1};
      varargin = varargin(2:end);
   end
   if isfield(g, 'parent'),
      if ~isfield(g.faces, 'centroids'),
         dispif(mrstVerbose, 'Adding (some) geometry to coarse grid.');
         g = computeCoarseGeometry(g);
      end
   elseif ~isfield(g.faces, 'centroids'),
      dispif(mrstVerbose, 'Adding geometry to grid.');
      g = computeGeometry(g);
   else
      error('');
   end
   v = ishold;hold on;
   coord = num2cell(g.faces.centroids(f,:), 1);
   h=text(coord{:}, num2str(f(:)), varargin{:});
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
      cg.faces.centroids(:,i) = accumarray(faceno, cg.parent.faces.centroids(cg.faces.fconn, i));
   end
   cg.faces.centroids = bsxfun(@rdivide, cg.faces.centroids, accumarray(faceno, cg.parent.faces.areas(cg.faces.fconn)));
end
