function state = upscaleState(coarsemodel, model, state)
    p = coarsemodel.G.partition;
    CG = coarsemodel.G;
    
    pvc = coarsemodel.operators.pv;
    pvf = model.operators.pv;
    
    counts = accumarray(p, 1);
    
    nph = size(state.s, 2);
    pvs = bsxfun(@times, state.s, pvf);
    
    % Calculate saturations based on new and old pore volume
    s = zeros(CG.cells.num, nph);
    for i = 1:nph
        s(:, i) = accumarray(p, pvs(:, i))./pvc;
    end
    state.s = s;
    
    % Average the pressure (not entirely correct for compressible systems,
    % but we won't start evaluating properties in here).
    state.pressure = accumarray(p, state.pressure)./counts;
    
    state.flux = zeros(CG.faces.num, nph);
end
