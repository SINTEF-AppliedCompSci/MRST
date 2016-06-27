function state = updateNextRS(model, state, problem, dx, drivingForces)
    
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
