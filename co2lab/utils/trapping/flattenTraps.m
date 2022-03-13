function Gt = flattenTraps(Gt, res)
% Given output res from trapAnalysis, produce a grid with flat top surfaces in trap areas.
%
% SYNOPSIS:
%   Gt = flattenTraps(Gt, res)
%
% PARAMETERS:
%   Gt  - Top surface grid
%
%   res - Output from trapAnalysis called on Gt.
%
% RETURNS:
%   Gt  - Gt where cells which are part of traps have had their z value
%   reduced to the level to the trap. Useful to make nice plots. This
%   function also adds the sortedCellNodes field to Gt if not present.
%
% SEE ALSO:
%   `interactiveTrapping`.

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
    for i = 1:max(res.traps)
        Gt.cells.z(res.traps == i) = res.trap_z(i);
    end
    z_spill_loc = zeros(Gt.cells.num, 1);
    
    ind = res.traps == 0;
    z_spill_loc(ind) = 0;
    z_spill_loc(~ind) = res.trap_z(res.traps((~ind)));

    cells     = find(z_spill_loc > 0); 
    eIX       = Gt.cells.facePos;      
    nn        = double(diff([Gt.cells.facePos(cells), ...
                             Gt.cells.facePos(cells + 1)], [], 2));
                         
    Gt.cells.sortedCellNodes = getSortedCellNodes(Gt);


    cn        = double(Gt.cells.sortedCellNodes(mcolon(eIX(cells), eIX(cells + 1) - 1), 1));
    zz        = rldecode(z_spill_loc(z_spill_loc>0),nn);
    Gt.nodes.z(cn) = zz;

    Gt = computeGeometry(Gt);
end
