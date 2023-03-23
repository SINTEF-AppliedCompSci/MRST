function grdecl = createPinchedColumn
%Create a single column containing a single pinched layer of thickess 0.01
%
% SYNOPSIS:
%   grdecl = createPinchedColumn
%
% PARAMETERS:
%   None.
%
% RETURNS:
%   grdecl - Corner-point specification suitable for passing to grid
%            constructor `processGRDECL`.
%
% SEE ALSO:
%   `processGRDECL`.

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


   cartDims = [1, 1, 3];

   [x, y, z] = ndgrid(0:1, 0:1, 0);
   coord     = repmat([x(:), y(:), z(:)], [1, 2]);

   d = 1.01;
   zcorn = ...
      [ 0 0 ; ...
        0 0 ; ...
        1 1 ; ...
        1 1 ; ...
              ...
        1 1 ; ...
        1 1 ; ...
        d d ; ...
        d d ; ...
              ...
        d d ; ...
        d d ; ...
        2 2 ; ...
        2 2 ; ...
      ];

   actnum = [ 1, 0, 1 ];
   pinch  = { 0.02, 'GAP', inf, 'TOPBOT', 'TOP' };

   grdecl = struct('cartDims', cartDims, ...
                   'PINCH'   , { pinch }, ...
                   'COORD'   , reshape(coord .', [], 1), ...
                   'ZCORN'   , reshape(zcorn .', [], 1), ...
                   'ACTNUM'  , reshape(actnum  , [], 1));
end
