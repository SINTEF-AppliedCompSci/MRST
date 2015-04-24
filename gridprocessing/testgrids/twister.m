function pt = twister(pt, varargin)
%Permutes x- and y-coordinates of nodes in a grid.
%
% SYNOPSIS:
%   pt = twister(pt)
%
% PARAMETERS:
%   pt - Coordinates to permute
%
% RETURNS:
%   pt - Modified coordinates.
%
% NOTE:
%   If 'pt' is the nodes.coords field of a 'grid_structure', then function
%   'twister' invalidates any pre-computed face areas and cell volues (and
%   other fields).  Consequently, function 'computeGeometry' must be called
%   *after* a call to function 'twister'.
%
% EXAMPLE:
%   G = cartGrid([30, 20]);
%   G.nodes.coords = twister(G.nodes.coords);
%   G = computeGeometry(G);
%   plotCellData(G, G.cells.volumes, 'EdgeColor', 'k'), colorbar
%
% SEE ALSO:
%   grid_structure, cartGrid.

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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

   if nargin == 1,
      value = 0.03;
   else
      value = varargin{1};
   end

   if isstruct(pt)       && isfield(pt,       'nodes') && ...
      isstruct(pt.nodes) && isfield(pt.nodes, 'coords'),
      % 'pt' is (presumably) a 'grid_structure'.  Handle by recursive call.
      pt.nodes.coords = twister(pt.nodes.coords, value);
   elseif isnumeric(pt) && any(size(pt,2) == [2,3]),
      % 'pt' is a coordinate array.  Map individual coordinates.

      % Coordinate mapping function.
      f  = @(x,y) value * sin(pi * x) .* sin(3 * (-pi/2 + pi*y));

      % Normalize coordinates to interval [0,1].
      xi = bsxfun(@rdivide,                    ...
                  bsxfun(@minus, pt, min(pt)), ...
                  max(pt) - min(pt));

      % Create new coordinates.
      pt(:,1:2) = bsxfun(@times, ...
                         [xi(:,1) + f(xi(:,1), xi(:,2)),  ...
                          xi(:,2) - f(xi(:,2), xi(:,1))], ...
                         max(pt(:,1:2)) - min(pt(:,1:2)));
   else
      error(msgid('DataType:Unknown'), ...
            'Don''t know how to handle input data type ''%s''', class(pt));
   end
end
