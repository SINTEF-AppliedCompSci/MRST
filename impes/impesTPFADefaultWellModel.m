function wdp = impesTPFADefaultWellModel(state, wells, fluid)
%Default well pressure model for impesTPFA solver
%
% SYNOPSIS:
%   wdp = impesTPFADefaultWellModel(state, wells, fluid)
%
% PARAMETERS:
%   state - Reservoir and well solution structure.
%
%   wells - Well data structure as defined by function 'addWell'.
%
%   fluid - Fluid data structure as defined by function 'initEclipseFluid'
%           or similar.
%
% RETURNS:
%   wdp - Well perforation gravity pressure adjusments.  One scalar value
%         for each perforation in each well.  Specifically, the pressure
%         difference (p_w + wdp - pc) may be inserted directly into
%         Peaceman's connection flux model.
%
% SEE ALSO:
%   impesTPFA.

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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


   if ~isempty(wells),
      g  = gravity;
      wc = vertcat(wells.cells);

      if norm(g) > 0,

         nperf = reshape(cellfun('prodofsize', { wells.cells }), [], 1);

         [rho, rho] = fluid.pvt(state.pressure(wc), state.z(wc, :));   %#ok

         rho = sparse(rldecode(1 : numel(nperf), nperf, 2), ...
                      1 : numel(wc), 1) * rho;
         rho = bsxfun(@rdivide, rho, nperf);
         rho = sum(rho .* vertcat(wells.compi), 2);

         wdp = norm(g) .* vertcat(wells.dZ) .* rldecode(rho, nperf);

      else

         wdp = zeros(size(wc));

      end

   else
      wdp = [];
   end
end
