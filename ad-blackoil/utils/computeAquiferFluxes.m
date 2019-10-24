function q = computeAquiferFluxes(model, state, dt)
    
    [p, sW] = model.getProps(state, 'pressure', 'water');
    p_aq = model.getProp(state, 'aquiferpressures');
    V_aq = model.getProp(state, 'aquifervolumes');

    aquifers  = model.aquifers;
    aquind  = model.aquind;
    
    alpha     = aquifers(:, aquind.alpha);
    J         = aquifers(:, aquind.J);
    conn      = aquifers(:, aquind.conn);
    depthconn = aquifers(:, aquind.depthconn);
    depthaq   = aquifers(:, aquind.depthaq);
    C         = aquifers(:, aquind.C);
    aquid     = aquifers(:, aquind.aquid);
    
    nconn = size(conn, 1);
    naq = max(aquid);
    aquid2conn = sparse(aquid, (1 : nconn)', 1, naq, nconn)';
    
    p_aq = aquid2conn*p_aq;
    V_aq = aquid2conn*V_aq;

    p  = p(conn);
    sW = sW(conn);
    
    fluid = model.fluid;
    pcOW = 0;
    ph = model.getPhaseNames();
    isw = (ph == 'W');
    if isfield(fluid, 'pcOW') && ~isempty(sW)
        pcOW = model.getProps(state, 'CapillaryPressure');
        pcOW = pcOW{isw}(conn);
    end

    b = model.getProps(state, 'ShrinkageFactors');
    bW = b{isw}(conn);
    rhoW = bW.*fluid.rhoWS;
    
    Tc = C.*V_aq./J;
    if dt == 0
        coef = 1;
    else
        coef = (1 - exp(-dt./Tc))./(dt./Tc);
    end
    
    g = model.gravity(3);
    q = alpha.*J.*(p_aq + pcOW - p + rhoW.*g.*(depthconn - depthaq)).*coef;
    
end
