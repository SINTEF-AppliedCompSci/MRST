function state = updateWellCellSaturationsExplicit(model, state, problem, dx, drivingForces)
    
    W = drivingForces.Wells;
    
    inj = vertcat(drivingForces.Wells.sign) > 0;
    
    active = inj;
    active = true(numel(W), 1);
%     perf2well = getPerforationToWellMapping(W);
%     active = inj(perf2well);
    
    cells = vertcat(W(active).cells);
    
    dt = problem.dt;
    
    b = state.bfactor(cells, :);
    b0 = state.bfactor0(cells, :);
    
    s = state.s(cells, :);
    s0 = state.s0(cells, :);

    
    intx = model.operators.internalConn;
    div = model.operators.C(:, cells)';
    
    
    pv = repmat(model.operators.pv(cells), 1, size(s, 2));
    
    cqs = vertcat(state.wellSol(active).cqs);
    
    
    divV = b.*(div*state.flux(intx, :));
    
    mass0 = s0.*b0;
    
    if model.disgas
        assert(model.water && model.gas && model.oil);
        GI = 3;
        OI = 2;

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
        
        rsnew = min(rsnew, rsSat);
%         state.rs(cells) = rsnew;
        if norm((rsnew - rs0)./rs0, inf) > 0.1;
            state.rs(cells) = rsnew;
        end
    end    
    state.s = cap(state.s);
    if 1
%         state.s(:, 1) = 1 - state.s(:, 2) - state.s(:, 3);
        ind = 1:size(state.s, 2);
        
        isfill = ind == 2;
        state.s(cells, isfill) = 1 - sum(state.s(cells, ~isfill), 2);
    else
        state.s = bsxfun(@rdivide, state.s, sum(state.s, 2));
    end
    
    above = state.s(:, 3) > 0;
    state.rs(above) = state.rsSat(above);
end
