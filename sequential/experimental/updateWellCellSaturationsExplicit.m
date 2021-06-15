function state = updateWellCellSaturationsExplicit(model, state, problem, dx, drivingForces)
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
    
    
    bdivV = b.*(div*state.flux(intx, :));
    
    mass0 = s0.*b0;
    
    if model.gas && model.disgas
        assert(model.water && model.gas && model.oil);
        GI = 3;
        OI = 2;

        rs = state.rs(cells);
        rs0 = state.rs0(cells);
        
        bdivV(:, GI) = bdivV(:, GI) + bdivV(:, OI).*rs;
        
        % mass(:, GI) = mass(:, GI) + mass(:, OI).*rs;
        mass0(:, GI) = mass0(:, GI) + mass0(:, OI).*rs0;
    end
    
    
    cap = @(x) min(max(x, 0), 1);
    
    dmass = dt*(cqs - bdivV)./pv;
    ds = dmass./b;
    
    mass = mass0 + dmass;
    
    
%     state.s(cells, :) = state.s(cells, :).*(b0./b) + ds;
    state.s(cells, :) = s0.*(b0./b) + ds;
    
%     s = state.s(cells, :);
    if model.gas && model.disgas
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
    
    if model.gas && model.disgas
        above = state.s(:, 3) > 0;
        state.rs(above) = state.rsSat(above);
    end
end
