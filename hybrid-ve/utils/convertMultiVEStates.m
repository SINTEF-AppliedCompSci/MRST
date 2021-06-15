function states = convertMultiVEStates(model, states_c)
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

    ns = numel(states_c);
    states = cell(ns, 1);
    for i = 1:ns
        states{i} = convertState(model, states_c{i});
        
        % add converted sGmax field (NB: not converted from
        % states_c{i}.sGmax, but rather simply obtained from higher
        % resolution sG).
        if i == 1
            states{i}.sGmax = states{i}.s(:,2);
        else
            states{i}.sGmax = max(states{i}.s(:,2), states{i-1}.sGmax);
        end
    
        % what about converting flux?
    end
end

function state = convertState(model, state_c)
    CG = model.G;
    p = CG.partition;
    
    t = CG.cells.topDepth(p);
    T = CG.parent.cells.topDepth;
    
    b = CG.cells.bottomDepth(p);
    B = CG.parent.cells.bottomDepth;
    
    h_c = state_c.s(:, 2).*CG.cells.height;
    h = h_c(p);
    
    sG = getGasSaturationFromHeight(T, t, B, b, h);
    
    state = state_c;
    g = norm(model.gravity);
    if isfield(state_c, 'rho')
        rhow = state_c.rho(p, 1);
        rhog = state_c.rho(p, 2);
    else
        if isprop(model, 'EOSModel')
            rhow = model.fluid.rhoOS;
        else
            rhow = model.fluid.rhoWS;
        end
        rhog = model.fluid.rhoGS;
        rhog = repmat(rhog, model.G.cells.num, 1);
        rhow = repmat(rhow, model.G.cells.num, 1);
    end
    isFine = CG.cells.discretization == 1;
    
    state_c.pressure(~isFine) = state_c.pressure(~isFine) + g.*rhow(~isFine).*CG.cells.height(~isFine)/2;
    cz = (T + B)/2;
    pressure = state_c.pressure(p);
    p_c = getWaterPressureFromHeight(cz, t, B, b, h, pressure, g, rhow(p), rhog(p));
    p_c(isFine(CG.partition)) = pressure(isFine(CG.partition));
    
    state.pressure = p_c;
    state.s = [1-sG, sG];
    
    if isprop(model, 'EOSModel')
        state.T = state.T(p);
        state.L = state.L(p);
        state.K = state.K(p, :);
        state.components = state_c.x(CG.partition, :).*state.s(:, 1) + state_c.y(CG.partition, :).*state.s(:, 2);
        pureLiquid = state.s(:, 1) == 1;
        pureVapor = state.s(:, 2) == 1;
        
        state.x = ~pureLiquid.*state_c.x(CG.partition, :) + pureLiquid.*state.components;
        state.y = ~pureVapor.*state_c.y(CG.partition, :) + pureVapor.*state.components;
    end
end
