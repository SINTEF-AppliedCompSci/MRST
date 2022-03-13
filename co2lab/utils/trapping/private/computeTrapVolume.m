function vol = computeTrapVolume(Gt, tstruct, poro, trap_ixs)
% Compute total pore volume of traps
%
% SYNOPSIS:
%   vol = computeTrapVolume(Gt, tstruct, poro, trap_ix)
%
% PARAMETERS:
%   Gt       - top surface grid
%   tstruct  - trap structure of 'Gt', as produced by 'trapAnalysis'
%   poro     - porosity (one value per cell in Gt)
%   trap_ixs - indices of traps for which volume should be computed
%
% RETURNS:
%   vol - total pore volumes of the requested traps
%
% SEE ALSO:
% `trapAnalysis`

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

  vol = zeros(numel(trap_ixs),1);
    
  for i = 1:numel(trap_ixs)
      % find index of all trap cells in this trap
      cix = find(tstruct.traps == trap_ixs(i));
      
      trap_heights = tstruct.trap_z(trap_ixs(i)) - Gt.cells.z(cix);
      trap_heights = min(trap_heights, Gt.cells.H(cix));
      
      cell_vols = trap_heights .* Gt.cells.volumes(cix) .* poro(cix);
      vol(i) = sum(cell_vols);
  end
end
