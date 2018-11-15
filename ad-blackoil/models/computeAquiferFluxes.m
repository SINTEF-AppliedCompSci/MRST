function q = computeAquiferFluxes(model, p, sW, state, dt)
    
    aquifers  = model.aquifers;

    alpha     = aquifers.alpha;
    J         = aquifers.J;    
    conn      = aquifers.conn;    
    depthconn = aquifers.depthconn;    
    depthaq   = aquifers.depthaq;    
    C         = aquifers.C;    
    aqid2conn = aquifers.aqui2conn;    
    
    p_aq = model.getProp(state, 'aquiferpressures');
    p_aq = aqid2conn*p_aq;

    V_aq = model.getProp(state, 'aquifervolumes');
    V_aq = aqid2conn*V_aq;

    p  = p(conn);
    sW = sW(conn);
    
    pcOW = 0;
    if isfield(fluid, 'pcOW') && ~isempty(sW)
        pcOW  = fluid.pcOW(sW, 'cellInx', conn);
    end

    bW     = fluid.bW(p_aq, 'cellInx', conn);
    rhoW   = bW.*fluid.rhoWS;
    
    T = C.*V_aq/J;
    coef = (1 - exp(-dt/Tc))./(dt/Tc);
    
    q = alpha.*J.*(p_aq + pcOW - p + rhoW*g*(depthconn - depthaq)).*coef;
    
end
