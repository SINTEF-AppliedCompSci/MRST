function checkBoundaryBasis(G, CG, CS, bc)
%Check structural integrity of coarse system versus boundary conditions.
%
% SYNOPSIS:
%   checkBoundaryBasis(G, CG, CS, bc)
%
% PARAMETERS:
%   G, CG - Grid and coarse grid, respectively.
%
%   CS    - Coarse-scale linear system structure as defined by function
%           'generateCoarseSystem'.
%
%   bc    - Boundary condition structure as defined by function 'addBC'.
%
% RETURNS:
%   Nothing.  However, a warning is printed if the set of boundary
%   conditions identified by 'bc' include external coarse boundaries for
%   which no basis functions have been computed.
%
% NOTE:
%   This is an internal function in the MsMFE implementation.  Its
%   existence, syntax and semantics may change without prior notice.
%
% SEE ALSO:
%   `solveIncompFlowMS`, `solveIncompFlowMSSpeedUp`.

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


   [nsub, sub] = subFaces(G, CG);

   ext = find(any(CG.faces.neighbors == 0, 2));
   cf  = ext(~ cellfun('isempty', CS.basis(ext)));

   p = cumsum([1; nsub]);
   f = sub(mcolon(p(cf), p(cf + 1) - 1));

   have_bdry_basis    = false([G.faces.num, 1]);
   have_bdry_basis(f) = true;

   if ~all(have_bdry_basis(bc.face)),
      warning(msgid('BC:NoBasis'), ...
             ['Boundary condition affects coarse interface for which ', ...
              'no flow basis function has been computed.\n'           , ...
              'Use option ''bc'' in ''generateCoarseSystem'' to '     , ...
              'compute flow basis functions for all interfaces '      , ...
              'influenced by boundary conditions.']);
   end
end
