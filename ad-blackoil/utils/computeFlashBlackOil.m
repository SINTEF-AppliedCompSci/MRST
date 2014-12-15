function state = computeFlashBlackOil(state, state0, model, status)
    % Update state based on vapoil/disgas flash
    fluid = model.fluid;
    
    disgas = model.disgas;
    vapoil = model.vapoil;
    
    p = model.getProp(state, 'p');
    p0 = model.getProp(state0, 'p');
    
    so = model.getProp(state, 'so');
    so0 = model.getProp(state0, 'so');
    
    sw = model.getProp(state, 'sw');
    % sw0 = model.getProp(state0, 'sw');
    
    sg = model.getProp(state, 'sg');
    sg0 = model.getProp(state0, 'sg');
    
    rs = model.getProp(state, 'rs');
    rs0 = model.getProp(state0, 'rs');
    
    rv = model.getProp(state, 'rv');
    rv0 = model.getProp(state0, 'rv');
    
    etol = sqrt(eps);
    % Detrmine staus of updated cells -----------------------------------------
    watOnly  = sw > 1-etol;

    % phase transitions sg <-> rs  --------------------------------------------
    if ~disgas
        rsSat0 = rs0;
        rsSat  = rsSat0; 
        gasPresent = true;
    else
        st1 = status{1};
        rsSat0 = fluid.rsSat(p0);
        rsSat  = fluid.rsSat(p);
        gasPresent = or(and( sg > 0, ~st1), watOnly); % Obvious case
        % Keep oil saturated if previous sg is sufficiently large:
        ix1 = and( sg < 0, sg0 > etol);
        gasPresent = or(gasPresent, ix1);
        % Set oil saturated if previous rs is sufficiently large
        ix2 = and( and(rs > rsSat*(1+etol), st1), rs0 > rsSat0*(1-etol) );
        assert(all(sg(ix2)==0))
        gasPresent = or(gasPresent, ix2);
    end
    ix = sg < 0;
    sw(ix) = sw(ix)./(1-sg(ix));
    so(ix) = so(ix)./(1-sg(ix));
    sg(ix) = 0;

    % phase transitions so <-> rv
    if ~vapoil
        oilPresent = true;
        rvSat0 = rv0;
        rvSat  = rvSat0;
    else
        st2 = status{2};
        rvSat0   = fluid.rvSat(p0);
        rvSat    = fluid.rvSat(p);
        oilPresent = or(and( so > 0, ~st2), watOnly); % Obvious case
        % Keep gas saturated if previous so is sufficiently large
        ix1 = and( so < 0, so0 > etol);
        oilPresent = or(oilPresent, ix1);
        % Set gas saturated if previous rv is sufficiently large
        ix2 = and( and(rv > rvSat*(1+etol), st2), rv0 > rvSat0*(1-etol) );
        assert(all(so(ix2)==0))
        oilPresent = or(oilPresent, ix2);
    end
    ix = so < 0;
    sw(ix) = sw(ix)./(1-so(ix));
    sg(ix) = sg(ix)./(1-so(ix));
    so(ix) = 0;

    % make sure sw >=0
    ix = sw < 0;
    so(ix) = so(ix)./(1-sw(ix));
    sg(ix) = sg(ix)./(1-sw(ix));
    sw(ix) = 0;

    % Update saturated r-values -----------------------------------------------
    rs(gasPresent) = rsSat(gasPresent);
    rv(oilPresent) = rvSat(oilPresent);

    % Update undersatured r-values
    rs(~gasPresent) = min(rsSat(~gasPresent), rs(~gasPresent));

    % Update state ------------------------------------------------------------
    if model.water
        state = model.setProp(state, 'sw', sw);
    end
    state = model.setProp(state, 'so', so);
    state = model.setProp(state, 'sg', sg);
    
    state = model.setProp(state, 'rs', max(rs, 0));
    state = model.setProp(state, 'rv', max(rv, 0));

    state.status = oilPresent + 2*gasPresent;
end