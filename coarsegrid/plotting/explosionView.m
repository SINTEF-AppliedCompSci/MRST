function explosionView(G, q, varargin)
%Plot a partition as an explosion view
%
% SYNOPSIS:
%   explosionView(G, p)
%   explosionView(G, p, r)
%   explosionView(G, p, r, o)
%   explosionView(G, p, r, o)
%   explosionView(G, p, r, o, a, 'pn', 'pv', ..)
%
% PARAMETERS:
%   G - Grid data structure.
%   p - partition vector (cells with zero value are not shown)
%   r - radius multiplier, will determine the relative distance each grid
%       block will be moved in the radial direction (r>0, default: 0.2)  
%   o - origin to use for the explosion view (default: grid center)
%   a - transparency value for faces, one value per coarse block 
%
% EXAMPLES:
%   1) Partition a 2D Cartesian grid
%   G = computeGeometry(cartGrid([10 10]));
%   p = partitionCartGrid(G.cartDims, [3 3]);
%   explosionView(G, p);
%
%   2) Partition a cup-shaped grid
%   x = linspace(-2,2,41);
%   G = tensorGrid(x,x,x);
%   G = computeGeometry(G);
%   c = G.cells.centroids;
%   r = c(:,1).^2 + c(:,2).^2+c(:,3).^2;
%   G = removeCells(G, (r>1) | (r<0.25) | (c(:,3)<0));
%   p = compressPartition(partitionUI(G,[5 5 4]));
%   explosionView(G, p, .4); view(3), axis tight

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

assert(any(G.griddim==[2 3]),'G must be a 2D or 3D grid');
assert(isfield(G.cells, 'centroids'), ...
   'Grid does not supply cell centroids. Please use "computeGeometry" first');

if nargin<3 || ~isnumeric(varargin{1})
   rmult = .2;
else
   rmult = varargin{1}; varargin=varargin(2:end);
   assert(numel(rmult)==1, 'radius multiplier must be a scalar');
end
if isempty(varargin) || ~isnumeric(varargin{1})
   o = mean(G.cells.centroids);
else
   o = varargin{1}; varargin=varargin(2:end);
   assert(numel(o)==G.griddim,'Origin must be a valid coordinate');
end
if isempty(varargin) || ~isnumeric(varargin{1})
   alpha = ones(max(q),1);
else
   alpha = varargin{1}; varargin=varargin(2:end);
   assert(numel(alpha)==max(q),'One transparency entry per block');
end
a = max(G.cells.centroids) - min(G.cells.centroids); a=a./max(a);
if G.griddim==2
   for i=1:max(q)
      g = extractSubgrid(G, q==i);
      c = mean(bsxfun(@minus, g.cells.centroids, o),1)./a;
      [th,r] = cart2pol(c(:,1),c(:,2));
      [x,y] = pol2cart(th, rmult*r);
      sgn = sign(mean(bsxfun(@minus,g.nodes.coords,o))) .* sign([x,y]);
      g.nodes.coords = ...
         bsxfun(@plus, g.nodes.coords, sgn.*[x,y]);
      plotCellData(g, i*ones(g.cells.num,1), 'EdgeColor', 'k', varargin{:});
   end
elseif G.griddim==3
   for i=1:max(q)
      g = extractSubgrid(G, q==i);
      c = mean(bsxfun(@minus, g.cells.centroids, o),1)./a;
      [th,phi,r] = cart2sph(c(:,1),c(:,2),c(:,3));
      [x,y,z] = sph2cart(th, phi, rmult*r);
      sgn = sign(mean(bsxfun(@minus,g.nodes.coords,o))) .* sign([x,y,z]);
      g.nodes.coords = ...
         bsxfun(@plus, g.nodes.coords, sgn.*[x,y,z].*a);
      bf = boundaryFaces(g);
      plotFaces(g, bf, i*ones(numel(bf),1), 'EdgeColor', 'k', ...
         'FaceAlpha', alpha(i), varargin{:});
    end
end

end