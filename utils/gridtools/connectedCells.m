function s = connectedCells(G, varargin)
%Compute connected components of grid cell subsets.
%
% SYNOPSIS:
%   s = connectedCells(G, c)
%
% PARAMETERS:
%   G - Grid data structure.
%
%   c - Cell subset.  List (array) of grid cell indices, or logical mask
%       into the cells of 'G'.
%
% RETURNS:
%   s - Connected component set.  An n-element `cell` array, one element for
%       each component, of grid cells.  Specifically, s{i} contains the
%       grid cell indices of component 'i'.
%
% SEE ALSO:
%   `grid_structure`, `dmperm`.

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


   c = (1 : G.cells.num) .';
   if nargin > 1 && ~isempty(varargin{1}),
      if isnumeric(varargin{1}),
         c = varargin{1};
      elseif ~islogical(varargin{1}) || numel(varargin{1}) ~= G.cells.num,
         error(['Cell subset must be array of indices or ', ...
                'logical mask into grid cells']);
      end
   end

   active        = false([G.cells.num, 1]);
   active(c)     = true;
   nactive       = double(sum(active));
   renum         = zeros([G.cells.num, 1]);
   renum(active) = 1 : nactive;
   renum         = [0; renum];

   % Create adjacency matrix.
   N = renum(G.faces.neighbors + 1) + 1;
   A = sparse([N(:,1); N(:,2); (1 : nactive+1).'], ...
              [N(:,2); N(:,1); (1 : nactive+1).'],  1 );
   A = A(2:end, 2:end); % Exclude inactive cells

   % Compute connected components.
   [p, r, r] = dmperm(A);                                              %#ok

   % Form output.  One CELL entry for each component.
   s = arrayfun(@(b,e) c(p(b : e)), r(1:end-1), r(2:end)-1, ...
                'UniformOutput', false);
end
