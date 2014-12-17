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
    
    if model.disgas
        rs = state.rs(cells);
        divV(:, 3) = divV(:, 3) + divV(:, 2)*rs;
    end
    
    
    cap = @(x) min(max(x, 0), 1);
    
    snew = cap((dt*(cqs - divV)./pv + s0.*b0)./b);
%     snew = cap(dt*(cqs - divV)./(pv.*b) + s0);
    if 0
        snew = bsxfun(@rdivide, snew, sum(snew, 2));
    else
        snew(:, 2) = 1 - snew(:, 1) - snew(:, 3);
    end
    
%     if s0(:, 3) > 0
%         return
%     end
    state.s(cells, :) = snew;
   
end
