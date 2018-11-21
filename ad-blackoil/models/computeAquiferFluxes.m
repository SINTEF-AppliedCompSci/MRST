function q = computeAquiferFluxes(model, p, sW, state, dt)
    
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
    
    p_aq = model.getProp(state, 'aquiferpressures');
    p_aq = aquid2conn*p_aq;

    V_aq = model.getProp(state, 'aquifervolumes');
    V_aq = aquid2conn*V_aq;

    p  = p(conn);
    sW = sW(conn);
    
    fluid = model.fluid;
    pcOW = 0;
    if isfield(fluid, 'pcOW') && ~isempty(sW)
        pcOW  = fluid.pcOW(sW, 'cellInx', conn);
    end

    bW     = fluid.bW(p_aq, 'cellInx', conn);
    rhoW   = bW.*fluid.rhoWS;
    
    Tc = C.*V_aq./J;
    coef = (1 - exp(-dt./Tc))./(dt./Tc);
    
    g = model.gravity(3);
    q = alpha.*J.*(p_aq + pcOW - p + rhoW.*g.*(depthconn - depthaq)).*coef;
    
end
