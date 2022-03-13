function bN = blockNeighbourship(N, p, varargin)
%Derive coarse-scale neighbourship from fine-scale information
%
% SYNOPSIS:
%   bN = blockNeighbourship(N, p)
%   bN = blockNeighbourship(N, p, f)
%
% PARAMETERS:
%   N  - Fine-scale neighbourship definition.  Often, but not necessarily
%        always, equal to the 'faces.neighbors' field of a grid structure.
%
%   p  - Partition vector defining coarse blocks.
%
%   f  - Optional boolean flag that, if set, includes the external boundary
%        into the inter-block connection set.  Default value: unset/false
%        whence external boundary is excluded from neighbourship relations.
%
% RETURNS:
%   bN - Coarse-scale (block) neighbourship definition (unique inter-block
%        connections).
%
% NOTE:
%   This function uses SORTROWS.
%
% SEE ALSO:
%   `grid_structure`, `generateCoarseGrid`, `sortrows`.

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


   if nargin > 2,
      assert (numel(varargin{1}) == 1, ...
              'Flag ''f'' must be scalar value convertible to LOGICAL.');
   end

   p  = [0; p];
   bN = p(N + 1);

   % Identify block connections.
   pick = bN(:,1) ~= bN(:,2);
   if (nargin == 2) || ~varargin{1},
      % Exclude boundary (i.e., pick internal connections exclusively).
      pick = pick & all(bN ~= 0, 2);
   end

   % Extract unique inter-block connections (symmetry handled elsewhere).
   bN = unique(sort(bN(pick, :), 2), 'rows');
end
