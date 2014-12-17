function state = updateWellCellSaturationsExplicit(model, state, problem, dx, drivingForces)
    
    W = drivingForces.Wells;
    
    inj = vertcat(drivingForces.Wells.sign) > 0;
%     perf2well = getPerforationToWellMapping(W);
%     active = inj(perf2well);
    
    cells = vertcat(W(inj).cells);
    
    dt = problem.dt;
    
    b = state.bfactor(cells, :);
    b0 = state.bfactor0(cells, :);
    
    s = state.s(cells, :);
    s0 = state.s0(cells, :);

    
    intx = model.operators.internalConn;
    div = model.operators.C(:, cells)';
    
    
    pv = repmat(model.operators.pv(cells), 1, size(s, 2));
    
    cqs = vertcat(state.wellSol(inj).cqs);
    
    
    divV = b.*(div*state.flux(intx, :));
    
    mass0 = s0.*b0;
    mass = s.*b;
    
    GI = 3;
    OI = 2;
    if model.disgas
        rs = state.rs(cells);
        rs0 = state.rs0(cells);
        
        divV(:, GI) = divV(:, GI) + divV(:, OI).*rs;
        
        % mass(:, GI) = mass(:, GI) + mass(:, OI).*rs;
        mass0(:, GI) = mass0(:, GI) + mass0(:, OI).*rs0;
    end
    
    
    cap = @(x) min(max(x, 0), 1);
    
    dmass = dt*(cqs - divV)./pv;
    ds = dmass./b;
    
    mass = mass0 + dmass;
    
    
    state.s(cells, :) = state.s(cells, :) + ds;
    
    if model.disgas
        rsSat = state.rsSat(cells);
        
        rsnew = mass(:, GI)./(b(:, OI).*s(:, OI));
        
        overflow = max(rsnew - rsSat, 0);
        state.s(cells, 3) = overflow.*(b(:, OI).*s(:, OI))./b(:, GI);
        
%         overflow
        
        rsnew = min(rsnew, rsSat);
        if norm(rsnew - state.rs(cells), inf) > 1;
%             disp([state.rs(cells), rsnew, rsSat])
            state.rs(cells) = rsnew;
        end
    end
%     [state.s(cells, :), sum(state.s(cells, :))]

%     ds
    
%     state.s(cells, :)
%     state.s(cells, :) = state.s(cells, :) + ds;
    
    state.s = cap(state.s);
%     state.s(:, 1) = 1 - state.s(:, 2) - state.s(:, 3);
    state.s(:, 2) = 1 - state.s(:, 1) - state.s(:, 3);
%     state.s = bsxfun(@rdivide, state.s, sum(state.s, 2));
    
    above = state.s(:, 3) > 0;
    state.rs(above) = state.rsSat(above);
    
    
    
%     snew = cap((dt*(cqs - divV)./pv + s0.*b0)./b);
%     snew = cap(dt*(cqs - divV)./(pv.*b) + s0);
%     if 0
%         snew = bsxfun(@rdivide, snew, sum(snew, 2));
%     else
%         snew(:, 2) = 1 - snew(:, 1) - snew(:, 3);
%     end
%     
%     if s0(:, 3) > 0
%         return
%     end
%     state.s(cells, :) = snew;
   
end
