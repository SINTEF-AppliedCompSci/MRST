function CGm = changeCoarseCenterToWellLocation(CGm,W)
% Designates coarse nodes to cells containing wells. For coarse blocks with
% multiple cells, it selects a random well location as coarse node

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
if CGm.parent.griddim>2
    CGm = setCentersByWells(CGm, W);
else
    wcells = [W.cells];
    for i = 1:numel(wcells)
        in = CGm.partition(wcells(i));
        CGm.cells.centers(in) = wcells(i);
    end
end
return

function CG = setCentersByWells(CG, W)
    for i = 1:numel(W)
        c = W(i).cells;
        local = unique(CG.partition(c));
        for j = 1:numel(local)
            cix = c(CG.partition(c) == local(j));
            mid = ceil(numel(cix)/2);
            CG.cells.centers(local(j)) = cix(mid);
        end
    end
return