function grdecl = pinchedNode()
%Define two-cell corner-point specification with single, pinched vertex
%
% SYNOPSIS:
%   grdecl = pinchedNode()
%
% PARAMETERS:
%   None.
%
% RETURNS:
%   grdecl - A GRDECL structure suitable for further processing by function
%            `processGRDECL`, refinement by function `refineGrdecl`, or
%            permanent storage on disk by function `writeGRDECL`.
%
% EXAMPLE:
%   grdecl = refineGrdecl(pinchedNode, [20, 40, 20]);
%   g      = computeGeometry(processGRDECL(grdecl));
%   fluid  = initSingleFluid('mu', 1, 'rho', 0);
%   rock   = struct('perm', ones([g.cells.num, 1]));
%   s      = computeMimeticIP(g, rock);
%   bc     = fluxside([], g, 'LEFT' ,  1);
%   bc     = fluxside(bc, g, 'RIGHT', -1);
%   x      = initState(g, [], 0);
%   x      = solveIncompFlow(x, g, s, fluid, 'bc', bc)
%   plotCellData(g, x.pressure, ...
%                'EdgeColor', 'k', 'EdgeAlpha', 0.15, 'FaceAlpha', 0.5)
%   view(3), colorbar
%
% SEE ALSO:
%   `processGRDECL`, `refineGrdecl`, `writeGRDECL`.

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


   coord = [         ...
      0 0 0 , 0 0 1; ...
      1 0 0 , 1 0 1; ...
      2 0 0 , 2 0 1; ...
                     ...
      0 1 0 , 0 1 1; ...
      1 1 0 , 1 1 1; ...
      2 1 0 , 2 1 1; ...
      ];

   zcorn = [     ...
      0 0   0 1; ...
      0 0   0 0; ...
                 ...
      1 1   1 1; ...
      1 1   1 1; ...
      ];

   grdecl = struct('cartDims', [2, 1, 1]               , ...
                   'COORD'   , reshape(coord .', [], 1), ...
                   'ZCORN'   , reshape(zcorn .', [], 1), ...
                   'ACTNUM'  , ones([2, 1], 'int32'));
end
