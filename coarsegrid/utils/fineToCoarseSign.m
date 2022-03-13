function sgn = fineToCoarseSign(cg)
%Compute sign change between fine faces and coarse faces.
%
% SYNOPSIS:
%   sgn = fineToCoarseSign(cg)
%
% PARAMETERS:
%   cg  - Coarse grid including parent (fine) grid. (New definition).
%
% RETURNS:
%   sgn - fine-to-coarse sign for fine faces in cg.faces.fconn.
%
% SEE ALSO:
%   `coarsenGeometry`, `coarsenBC`, `coarsenFlux`.

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


   assert (isfield(cg, 'parent'), ...
          ['Field ''parent''missing in coarse grid.  ', ...
           'Did you really supply a coarse grid?']);

   g      = cg.parent;

   % Find sign relation between fine faces and coarse faces: the sign of
   % the flux is positive if there is a flux from g.faces.neighbors(:,1) to
   % g.faces.neighbors(:,2).  The contribution to the coarse flux is
   % positive if it is from the coarse block cg.faces.neighbors(:,1) to
   % cg.faces.neighbors(:,1).

   faceno = rldecode(1:cg.faces.num, diff(cg.faces.connPos), 2)';
   p      = [0; cg.partition];

   % First fine cell
   c1     = g.faces.neighbors(cg.faces.fconn,1);

   % First coarse block
   b1     = cg.faces.neighbors(faceno,1);

   sgn    = 2*( p(c1+1)==b1 ) -1;
end
