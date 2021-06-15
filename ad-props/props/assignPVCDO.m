function f = assignPVCDO(f, pvcdo, reg)
    [f.bO, f.muO] = getFunctions(pvcdo, reg);
end

function [bO, muO] = getFunctions(PVCDO, reg)
    [bO, muO] = deal(cell(1, reg.pvt));
    
    for i = 1:reg.pvt
        pvcdo = PVCDO(i, :);
        
        por  = pvcdo(1); % ref pres
        bor  = pvcdo(2); % ref fvf
        co   = pvcdo(3); % compress
        muor = pvcdo(4); % ref visc
        vbo  = pvcdo(5); % viscosibility

        bO{i}  = @(po) exp(co.*(po-por))./bor;
        muO{i} = @(po) muor.*exp(vbo.*(po-por));
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


