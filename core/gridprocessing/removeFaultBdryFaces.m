function [flt, keep] = removeFaultBdryFaces(flt, G)
%Remove fault faces on boundary
%
% SYNOPSIS:
%    flt       = removeFaultBdryFaces(flt, G)
%   [flt, msk] = removeFaultBdryFaces(flt, G)
%
% PARAMETERS:
%   flt - A fault structure as defined by function `processFaults`.
%
%   G   - A grid structure representing a discretised reservoir geometry.
%
% RETURNS:
%   flt - An updated fault structure.  Faults that no longer have any
%         constituent faces from the grid are removed.
%
%   msk - Logical mask into original `flt` that indicates whether or not
%         the corresponding fault structure was still active after removing
%         its boundary faces.  Specifically, `msk(i)==true` if fault `i`
%         still has active faces after removing boundary contributions.
%         Moreover, its position in the updated `flt` structure is `pos(i)`
%         where `pos=cumsum(double(msk))`.
%
% SEE ALSO:
%   `processFaults`

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


   nflt = numel(flt);
   ff   = vertcat(flt.faces);
   nff  = vertcat(flt.numf);

   % Count number of boundary faces, i.e. those faces for which one
   % neighbour cell is zero, for each fault declared in 'flt'.
   %
   isbdry = any(G.faces.neighbors(ff, :) == 0, 2);
   nbdry  = accumarray(rldecode(1 : nflt, nff, 2) .', isbdry, [nflt, 1]);

   % Position vector into
   p = cumsum([1; nff - nbdry]);

   % Keep those faults that have any faces after removing boundary faces.
   keep = diff(p) > 0;

   % Partition the remaining fault faces amongst the active faults.
   tmp = ff(~isbdry);
   tmp = arrayfun(@(i) tmp(p(i) : p(i+1) - 1), ...
                  find(keep), 'UniformOutput', false);

   [flt(keep).faces] = tmp{:};

   % Count number of faces of each remaining fault.
   tmp = cellfun(@numel, tmp, 'UniformOutput', false);
   [flt(keep).numf] = tmp{:};

   % Exclude empty faults (i.e., those faults that are no longer
   % represented by any faces from the grid structure).
   flt(~keep) = [];
end
