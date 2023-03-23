function grdecl = pinchMiddleCell(k)
%Create corner-point descriptions with variable number of pinched nodes
%
% SYNOPSIS:
%   grdecl = pinchMiddleCell(n)
%
% PARAMETERS:
%   n      - Number of middle cell pinched corners.  Must be between zero
%            and four, inclusive.  OPTIONAL.  Default value: n = 1.
%
% RETURNS:
%   grdecl - Cell array of corner-point descriptions suitable for
%            subsequent processing using function `processGRDECL`.  Array
%            size is `nchoosek(4,n)`, one element for each combination of
%            selecting `n` elements from four possibilities.
%
% NOTE:
%   This is an extension of the test case defined by function `pinchedNode`
%
% SEE ALSO:
%   `pinchedNode`, `processGRDECL`, `nchoosek`.

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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


   [x, y, z] = ndgrid(0 : 1, 0 : 1, 0);
   coord     = repmat([x(:), y(:), z(:)], [1, 2]);

   if nargin < 1, k = 1; end

   assert (isnumeric(k) && (numel(k) == 1) && (0 <= k) && (k <= 4)  , ...
          ['Input parameter ''n'' (number of pinched corners) must ', ...
           'be a single integer between zero and four inclusive.']);

   v = nchoosek(1 : 4, k);

   grdecl = cell([size(v,1), 1]);

   for r = 1 : size(v, 1),
      i = v(r,:);

      grdecl{r} = struct('COORD'   , reshape(coord     .', [], 1), ...
                         'ZCORN'   , reshape(corners(i).', [], 1), ...
                         'cartDims', [1, 1, 3]);
   end
end

%--------------------------------------------------------------------------

function zcorn = corners(i)
   cellbound    = ones([1, 4]);
   cellbound(i) = 2;

   zcorn = [       ...
      0  0  0  0,  ...
      cellbound ,  ...
                   ...
      cellbound ,  ...
      2  2  2  2,  ...
                   ...
      2  2  2  2,  ...
      3  3  3  3,  ...
      ];
end
