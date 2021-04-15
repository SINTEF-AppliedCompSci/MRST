function varargout = outlineCoarseGrid(G, p, varargin)
%Impose outline of coarse grid on existing grid plot.
%
% SYNOPSIS:
%       outlineCoarseGrid(G, p)
%       outlineCoarseGrid(G, p, 'pn1', pv1, ...)
%       outlineCoarseGrid(G, p, c, 'pn1', 'pv1', ...)
%   h = outlineCoarseGrid(...)
%
% PARAMETERS:
%   G       - Grid data structure.
%
%   p       - Coarse grid partition vector as defined by (e.g) partitionUI
%             and processPartition.
%
%   c       - color, works only if G.griddim==2
%
% KEYWORD ARGUMENTS:
%
%   'Any'   - Additional keyword arguments will be passed directly on to
%             function `patch` meaning all properties supported by `patch`
%             are valid.
%
% RETURNS:
%   h - Handle to polygonal patch structure as defined by function
%       plotFaces. Only returned if requested.
%
% EXAMPLE:
%   G  = cartGrid([8, 8, 2]);
%   p  = partitionUI(G, [2, 2, 1]);
%   % plot fine grid:
%   plotGrid(G, 'faceColor', 'none'); view(3);
%   % outline coarse grid on fine grid:
%   outlineCoarseGrid(G, p);
%
% SEE ALSO:
%   `plotFaces`, `plotGrid`, `patch`

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


   assert (sum(G.griddim == [2, 3]) == 1);

   if G.griddim == 2,
      % Include outer boundary in 2D.
      mask = true([G.faces.num, 1]);
   else
      assert(~mod(nargin,2), ...
        'Function must be called with an even number of arguments');
      mask = all(G.faces.neighbors > 0, 2);
   end

   p = [0; p];
   N = p(G.faces.neighbors + 1);
   N = find((N(:,1) ~= N(:,2)) & mask);

   colour = [0.85, 0.15, 0.85];
   if G.griddim == 2,
      if (nargin > 2) && mod(nargin,2),
         colour = varargin{1};
         varargin = varargin(2:end);
      end
      h = plotEdges(G, N, colour, varargin{:});
   else
      h = plotFaces(G, N, colour, varargin{:});
   end

   if nargout > 0, varargout{1} = h; end
end

%--------------------------------------------------------------------------

function h = plotEdges(G, e, c, varargin)
   assert (all(diff([G.faces.nodePos(  e  ), ...
                     G.faces.nodePos(e + 1)], [], 2) == 2));

   fn = G.faces.nodes(mcolon(G.faces.nodePos(  e  ), ...
                             G.faces.nodePos(e + 1) - 1));
   fc = G.nodes.coords(fn,:);

   xd = reshape(fc(:,1), 2, []);
   yd = reshape(fc(:,2), 2, []);
   h  = line(xd, yd, 'Color', c, 'LineWidth', 2, varargin{:});
end
