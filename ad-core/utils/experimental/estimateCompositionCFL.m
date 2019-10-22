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
    q = value(model.getProp(state,'ComponentTotalFlux')');
    compMass = value(model.getProp(state, 'ComponentTotalMass')');
    totMass = sum(compMass, 2);
    % Ignore small masses 
    bad = compMass./totMass < 1e-7;
    compMass(bad) = 1;
    massFace = upstream(model, compMass, flag, xflow, l, r);
    
    rate_face = abs(q)./massFace;    
    % Accumulate into cell if flow is outgoing, or we have any kind of
    % cross-flow.
    nc = model.G.cells.num;
    ncomp = model.getNumberOfComponents();
    cfl = zeros(nc, ncomp);
    
    for i = 1:ncomp
        rate_cell = accumarray(l, rate_face(:, i).*( flag | xflow), [nc, 1]) +...
                    accumarray(r, rate_face(:, i).*(~flag | xflow), [nc, 1]);
        cfl(:, i) = dt.*rate_cell;
    end
    % Guard against cells without inflow
    cfl(~isfinite(cfl)) = 0;
end

function df_face = upstream(model, F, flag, xflow, l, r)
    df_face = model.operators.faceUpstr(flag, F);
    df_face(xflow) = max(abs(F(l(xflow))), abs(F(r(xflow))));
end
