function f = assignPVTG(f, pvtg, reg)
    [f.bG, f.muG, f.rvSat] = getFunctions(pvtg, reg);
end

function [bG, muG, rvSat] = getFunctions(PVTG, reg)
    [bG, muG, rvSat] = deal(cell(1, reg.pvt));
    
    for i = 1:reg.pvt
        pvtg = PVTG{i};
        
        rv = pvtg.data(pvtg.pos(1:end-1),1);
        p_vap = pvtg.key;
        
        % Cap endpoint (rv -> 0 at surface conditions)
        p_t = [0; p_vap];
        rv_t = [0; rv];

        bg = pvtg;
        bg.data = [bg.data(:,1), 1./bg.data(:,2)];
        bg = preprocessTablePVT(bg);
        
        mug = pvtg;
        mug.data = [mug.data(:,1), mug.data(:,3)];
        mug = preprocessTablePVT(mug);
        
        m = reg.pvtMethodGas;
        bG{i} = @(pg, rv, flag) interpPVT(bg, rv, pg, flag, m, reg.useMex);
        muG{i} = @(pg, rv, flag) interpPVT(mug, rv, pg, flag, m, reg.useMex);
        rvSat{i} = @(pg) reg.interp1d(p_t, rv_t, pg);
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
