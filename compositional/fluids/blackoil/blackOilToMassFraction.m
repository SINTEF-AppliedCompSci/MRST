function [xo, xg, yo, yg, zo, zg, rhoO, rhoG, muO, muG, freeGas] = blackOilToMassFraction(model, p, sO, sG, rs, rv)
% Evaluate properties

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

    fluid = model.fluid;
    rhoGS = fluid.rhoGS;
    rhoOS = fluid.rhoOS;
    
    disgas = model.disgas;
    if disgas
        isSatRs = rs >= fluid.rsSat(p);
        bO = fluid.bO(p, rs, isSatRs);

        muO = fluid.muO(p, rs, isSatRs);
        rhoO = bO.*(rs*rhoGS + rhoOS);
    else
        bO = fluid.bO(p);
        muO = fluid.muO(p);
        rhoO = bO.*rhoOS;
    end

    bG = fluid.bG(p);
    muG = fluid.muG(p);
    
    rhoG = bG.*rhoGS;
    zg = (rhoG.*sG + disgas.*rhoGS.*rs.*bO.*sO)./(rhoO.*sO + rhoG.*sG);
    zo = 1 - zg;
    
    xo = bO.*rhoOS./rhoO;
    xg = disgas.*bO.*rhoGS.*rs./rhoO;
    yg = 1;
    yo = 0;
    
    if model.disgas
        freeGas = isSatRs;
    else
        freeGas = zg > 0;
    end
end
