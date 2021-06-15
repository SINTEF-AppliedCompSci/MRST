function f = assignPVTO(f, pvto, reg)
    [f.bO, f.muO, f.rsSat, f.pb] = getFunctions(pvto, reg);
end

function [bO, muO, rsSat, pb] = getFunctions(PVTO, reg)
    [bO, muO, rsSat, pb] = deal(cell(1, reg.pvt));
    
    for i = 1:reg.pvt
        pvto = PVTO{i};
        % Transform B into reciprocal b = 1/b
        bo = pvto;
        bo.data = [bo.data(:,1), 1./bo.data(:,2)];
        bo = preprocessTablePVT(bo);

        % Viscosity
        muo = pvto;
        muo.data = [muo.data(:,1), muo.data(:,3)];
        muo = preprocessTablePVT(muo);
        % Bubble point vs rs table
        p_bub = pvto.data(pvto.pos(1:end-1),1);
        rs = pvto.key;
        % Cap endpoint to avoid interpolating into negative rs for low
        % pressures. The Rs should by definition go to zero as the pressure
        % goes towards surface conditions, and by extension, zero.
        p_t = [0; p_bub];
        rs_t = [0; rs];
        
        m = reg.pvtMethodOil;
        bO{i} = @(po, rs, flag) interpPVT(bo, po, rs, flag, m, reg.useMex);
        muO{i} = @(po, rs, flag) interpPVT(muo, po, rs, flag, m, reg.useMex);
        rsSat{i} = @(po) reg.interp1d(p_t, rs_t, po);
        pb{i} = @(rsi) reg.interp1d([-1; rs], [p_bub(1); p_bub], rsi);
    end
end

% 
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
