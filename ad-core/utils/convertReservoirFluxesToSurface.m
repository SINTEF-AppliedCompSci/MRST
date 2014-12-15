function states = convertReservoirFluxesToSurface(model, states)
    if numel(states) > 1
        assert(iscell(states));
    else
        states = {states};
    end
    
    if ~isfield(states{1}, 'bfactor')
        error(['Missing bfactors in state - were the states produced from',...
               ' a model with the ''''extraStateOutput'''' flag enabled']);
    end
    
    upstr = model.operators.faceUpstr;
    intx  = model.operators.internalConn;
    
    for i = 1:numel(states)
        state = states{i};
        
        [m, n] = size(state.flux);
        state.surfaceFlux = zeros(m, n);
        for j = 1:n
            up = state.upstreamFlag(:, j);
            b = state.bfactor(:, j);
            v = state.flux(intx, j);
            
            state.surfaceFlux(intx, j) = upstr(up, b).*v;
        end
        states{i} = state;
    end
end