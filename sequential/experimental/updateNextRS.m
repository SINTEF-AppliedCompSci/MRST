function state = updateNextRS(model, state, problem, dx, drivingForces)
%Undocumented Utility Function

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
    
    W = drivingForces.Wells;
    wc = vertcat(W.cells);
    
    dt = problem.dt;
    
    b = state.bfactor;
    
    s = state.s;

    
    intx = model.operators.internalConn;
    
    pv = model.operators.pv;
    
    cqs = vertcat(state.wellSol.cqs);
    
    
    
    
    GI = 3;
    OI = 2;
    
    qG = cqs(:, GI);
    upO = state.upstreamFlag(:, OI);
    
    bO = b(:, OI);
    sO = s(:, OI);
    
    
    rs0 = state.rs0;
    rs = state.rs;
    % SIGN HERE???
    gas = -model.operators.Div(model.operators.faceUpstr(upO, rs0.*bO).*state.flux(intx, OI));
    gas(wc) = gas(wc) + qG;
    
    state.rs = rs0 + dt*gas./(bO.*max(sO, 0.001).*pv);
    
%     state.rs(1)./rs0(1)
%     state.rs = rs + dt*(qG - RsbOvO)./(bO.*max(sO, 0.001).*pv);
    
    isSat = s(:, GI) > 0;
    

    
    state.rs(isSat) = model.fluid.rsSat(state.pressure(isSat));
    
    
%     if model.gas && model.disgas
%         assert(model.water && model.gas && model.oil);
%         GI = 3;
%         OI = 2;
% 
%         rs = state.rs(cells);
%         rs0 = state.rs0(cells);
%         
%         bdivV(:, GI) = bdivV(:, GI) + bdivV(:, OI).*rs;
%         
%         % mass(:, GI) = mass(:, GI) + mass(:, OI).*rs;
%         mass0(:, GI) = mass0(:, GI) + mass0(:, OI).*rs0;
%     end
    
%     dmass = dt*(cqs - bdivV)./pv;

end
