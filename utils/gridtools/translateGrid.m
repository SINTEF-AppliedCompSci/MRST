function G = translateGrid(G, dx)
%Move all grid coordinates according to particular translation
%
% SYNOPSIS:
%   G = translateGrid(G, dx)
%
% PARAMETERS:
%   G  - MRST grid as outlined in `grid_structure`.
%
%   dx - User-specified constant translation vector.  Must be a single
%        vector that will applied to all grid coordinates (i.e., nodes and
%        centroids).  Number of vector components must be equal to the
%        number of vertex coordinates (i.e., `size(G.nodes.coords, 2)`).
%
% RETURNS:
%   G - Updated grid structure.

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

   assert (numel(dx) == size(G.nodes.coords, 2), ...
          ['Translation vector must specify one component ', ...
           'in each coordinate direction.\nExpected %d components, ', ...
           'but got %d instead.'], size(G.nodes.coords, 2), numel(dx));

   dx        = reshape(dx, 1, []);  % BSXFUN needs row vector.
   translate = @(v) bsxfun(@plus, v, dx);

   G.nodes.coords = translate(G.nodes.coords);

   for s = { 'cells', 'faces' },
      if isfield(G.(s{1}), 'centroids'),
         G.(s{1}).centroids = translate(G.(s{1}).centroids);
      end
   end

   G.type = [ G.type, { mfilename } ];
end
