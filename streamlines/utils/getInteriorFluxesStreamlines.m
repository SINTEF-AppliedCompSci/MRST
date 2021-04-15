function flux = getInteriorFluxesStreamlines(G, state, pvol, reverse)
%Undocumented Utility Function

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

   if nargin < 4
       reverse = false;
   end
   if size(state.flux, 2) > 1
      state.flux = sum(state.flux, 2);
   end

   if reverse
      state.flux = -state.flux;
   end
   % Make array face fluxes for each cell in grid (Not outer).
   cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
   cf     = G.cells.faces;
   flux   = accumarray([cellNo, cf(:,2)], state.flux(cf(:,1)));
   flux   = bsxfun(@rdivide, flux, pvol);
end
