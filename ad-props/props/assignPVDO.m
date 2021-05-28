function f = assignPVDO(f, pvdo, reg)
    [f.bO, f.muO] = getFunctions(pvdo, reg);
end


function [bO, muO] = getFunctions(PVDO, reg)
    [bO, muO] = deal(cell(1, reg.pvt));
    
    for i = 1:reg.pvt
        pvdo = PVDO{i};
        p = pvdo(:, 1);
        BO = pvdo(:, 2);
        muo = pvdo(:, 3);
        bO{i}  = @(po) reg.interp1d(p, 1./BO, po);
        muO{i} = @(po) reg.interp1d(p, muo, po);
    end
end

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
