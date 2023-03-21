function grdecl = raisedColumn
%Create corner-point description of 2-by-1-by-2 grid with one fault
%
% SYNOPSIS:
%   grdecl = raisedColum
%
% PARAMETERS:
%   None.
%
% RETURNS:
%   grdecl - Corner-point description.
%
% NOTE:
%   This is primarily intended as a test grid in the development of an
%   edge-conforming corner-point processing routine.
%
% SEE ALSO:
%   `simpleGrdecl`, `processGRDECL`.

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


   [x, y, z] = ndgrid(0 : 2, 0 : 1, 0);
   coord     = repmat([x(:), y(:), z(:)], [1, 2]);

   zcorn = ...
      [  0 0 , 0   0   ; ...
         0 0 , 0   0   ; ...
         1 1 , 0.5 0.5 ; ...
         1 1 , 0.5 0.5 ; ...
                         ...
         1 1 , 0.5 0.5 ; ...
         1 1 , 0.5 0.5 ; ...
         2 2 , 1.5 1.5 ; ...
         2 2 , 1.5 1.5 ; ...
      ];

   grdecl = struct('cartDims', [2, 1, 2], ...
                   'COORD'   , reshape(coord .', [], 1), ...
                   'ZCORN'   , reshape(zcorn .', [], 1));
end
