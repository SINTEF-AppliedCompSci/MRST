function cfl = estimateCompositionCFL(model, state, dt, varargin)
    opt = struct('forces', []);
    opt = merge_options(opt, varargin{:});
    
    pv = model.operators.pv;
    internal = model.operators.internalConn;
    v = state.flux(internal, :);
    vT = sum(v, 2);
    
    xflow = ~(all(v >= 0, 2) | all(v <= 0, 2));
    
    l = model.operators.N(:, 1);
    r = model.operators.N(:, 2);
    
    flag = vT > 0;
    
    oix = 1+model.water;
    gix = oix+1;
    
    q = state.componentFluxes(internal, :);
    so = state.s(:, oix);
    sg = state.s(:, gix);
    
    rhoo = state.rho(:, oix);
    rhog = state.rho(:, gix);
    
    X = model.EOSModel.getMassFraction(state.x);
    Y = model.EOSModel.getMassFraction(state.y);
    
    compMass = X.*rhoo.*so + Y.*rhog.*sg;
    totMass = rhoo.*so + rhog.*sg;
    % Ignore small masses 
    bad = compMass./totMass < 1e-7;
    compMass(bad) = 1;
    compMass(so + sg == 0, :) = 0;
    massFace = upstream(model, compMass, flag, xflow, l, r);
    
    rate_face = abs(q)./massFace;    
    % Accumulate into cell if flow is outgoing, or we have any kind of
    % cross-flow.
    [nc, ncomp] = size(state.components);
    cfl = zeros(nc, ncomp);
    
    for i = 1:ncomp
        rate_cell = accumarray(l, rate_face(:, i).*( flag | xflow), [nc, 1]) +...
                    accumarray(r, rate_face(:, i).*(~flag | xflow), [nc, 1]);
        cfl(:, i) = (dt./pv).*rate_cell;
    end
    % Guard against cells without inflow
    cfl(~isfinite(cfl)) = 0;
end

function df_face = upstream(model, F, flag, xflow, l, r)
    df_face = model.operators.faceUpstr(flag, F);
    df_face(xflow) = max(abs(F(l(xflow))), abs(F(r(xflow))));
end
