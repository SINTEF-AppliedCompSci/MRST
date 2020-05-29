function f = assignPVTO(f, pvto, reg)
    [f.bO, f.muO, f.rsSat] = getFunctions(pvto, reg);
end

function [bO, muO, rsSat] = getFunctions(PVTO, reg)
    [bO, muO, rsSat] = deal(cell(1, reg.pvt));
    
    for i = 1:reg.pvt
        pvto = PVTO{i};
        
        p_bub = pvto.data(pvto.pos(1:end-1),1);
        rs = pvto.key;
        
        bo = pvto;
        bo.data = [bo.data(:,1), 1./bo.data(:,2)];
        
        muo = pvto;
        muo.data = [muo.data(:,1), muo.data(:,3)];
        
        bO{i} = @(po, rs, flag) interpPVT(bo, po, rs, flag);
        muO{i} = @(po, rs, flag) interpPVT(muo, po, rs, flag);
        rsSat{i} = @(po) reg.interp1d(p_bub, rs, po);
    end
end

% 
%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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
