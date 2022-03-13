function sc = coarse_sat(s, p, pv, nc)
% Converts a fine saturation field to a coarse saturation field, weighted
% with pore volumes.
%
% SYNOPSIS:
%   sc = coarse_sat(s, partition, pv, nc)
%
% DESCRIPTION:
%   Converts a fine saturation field to a coarse saturation field, weighted
%   with pore volumes.
%
% PARAMETERS:
%   s           - Fine saturation field
%   partitition - A partition vector representing the coarse grid
%   pv          - the pore volumes of the fine grid cells.
%                 pv=poreVolume(G, rock) at the fine scale.
%   nc          - Number of coarse grid blocks.
%
% RETURNS:
%   sc - Coarse saturation field, length nc.

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


   sc = accumarray(p, s(:,1) .* pv, [nc, 1]) ./ ...
        accumarray(p,           pv, [nc, 1]);
end
