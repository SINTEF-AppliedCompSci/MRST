function CS = assignBasisFuncs(CS, V, P)
%Update reservoir coarse system basis functions following 'evalBasisFunc'.
%
% SYNOPSIS:
%   CS = assignBasisFuncs(CS, V, P)
%
% PARAMETERS:
%   CS   - Original coarse system.  Must have allocated fields '.basis'
%          and '.basisP' of correct sizes.
%
%   V, P - New basis function values for flux (V) and pressure (P).
%          Assumed to be the return values from function 'evalBasisFunc'.
%
% RETURNS:
%   CS - Updated reservoir coarse system structure with new values for flux
%        and pressure basis functions.
%
% SEE ALSO:
%   `generateCoarseSystem`, `evalBasisFunc`.

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


   % Extract coarse faces for which new basis function values are defined.
   f = cellfun(@(x) x{3}, V);

   CS.basis (f) = V;
   CS.basisP(f) = P;
end
