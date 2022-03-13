function rock_cg = convertRock2Coarse(G, CG, rock)
%Create coarse-scale porosity field for solving transport equation.
%
% SYNOPSIS:
%   rock_cg = convertRock2Coarse(G, CG, rock)
%
% DESCRIPTION:
%   Computes coarse scale porosity as a volume-weighted average of
%   fine scale porosites within each coarse block.
%
%   Does not produce a coarse scale equivalent of rock.perm, whence the
%   coarse scale rock data is unsuited to models involving gravity forces.
%
% REQUIRED PARAMETERS:
%   G    - A grid_structure.
%
%   CG   - Coarse grid data sturcture as defined by function
%         'generateCoarseGrid'.
%
%   rock - Rock data structure with an accumulated porosity field.
%
% RETURNS:
%   rock_cg - Rock data structure defined on coarse grid, 'CG'.

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


   rock_cg.poro = accumarray(CG.partition, poreVolume(G, rock)) ./ ...
                  accumarray(CG.partition, G.cells.volumes);
end
