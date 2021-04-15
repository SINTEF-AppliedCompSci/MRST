function cfl = estimateCompositionCFL(model, state, dt, varargin)
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

    opt = struct('forces', [], 'useInflow', false);
    opt = merge_options(opt, varargin{:});
    
    internal = model.operators.internalConn;
    v = state.flux(internal, :);
    vT = sum(v, 2);
    
    xflow = ~(all(v >= 0, 2) | all(v <= 0, 2));
    
    l = model.operators.N(:, 1);
    r = model.operators.N(:, 2);
    
    flag = vT > 0;
    q = value(model.getProp(state,'ComponentTotalFlux')');
    compMass = value(model.getProp(state, 'ComponentTotalMass')');
    totMass = sum(compMass, 2);
    % Ignore small masses 
    bad = bsxfun(@rdivide, compMass, totMass) < 1e-7;
    compMass(bad) = 1;
    massFace = upstream(model, compMass, flag, xflow, l, r);
    
    rate_face = abs(q)./massFace;    
    % Accumulate into cell if flow is outgoing, or we have any kind of
    % cross-flow.
    nc = model.G.cells.num;
    ncomp = model.getNumberOfComponents();
    cfl = zeros(nc, ncomp);
    
    for i = 1:ncomp
        if opt.useInflow
            % Inflow faces is used to estimate throughput
            f = ~flag;
        else
            % Outflow faces is used to estimate throughput
            f = flag;
        end
        rate_cell = accumarray(l, rate_face(:, i).*( f | xflow), [nc, 1]) +...
                    accumarray(r, rate_face(:, i).*(~f | xflow), [nc, 1]);
        cfl(:, i) = dt.*rate_cell;
    end
    % Guard against cells without inflow
    cfl(~isfinite(cfl)) = 0;
end

function df_face = upstream(model, F, flag, xflow, l, r)
    df_face = model.operators.faceUpstr(flag, F);
    df_face(xflow) = max(abs(F(l(xflow))), abs(F(r(xflow))));
end
