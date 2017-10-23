function [problem, state] = equationsTracer(state0, state, model, dt, drivingForces, varargin)

    nt = model.getNumberOfTracers();
    
    [tracers, tracers0] = deal(cell(nt, 1));
    primaryVars = model.tracerNames;
    [tracers{:}] = model.getProps(state, primaryVars{:}, 'wellSol');
    [tracers0{:}] = model.getProps(state0, primaryVars{:});
    
    [tracers{:}] = initVariablesADI(tracers{:});
    
    [vT, qT] = getFlux(model, state); 
    
    [eqs, types] = deal(cell(1, nt));
    names = primaryVars;
    op = model.operators;
    flag = vT > 0;
    
    W = drivingForces.W;
    wc = vertcat(W.cells);
    p2w = getPerforationToWellMapping(W);
    comp = vertcat(W.tracer);
    comp = comp(p2w, :);
    
    inj = qT > 0;
    for i = 1:nt
        t = tracers{i};
        t0 = tracers0{i};
        vi = op.faceUpstr(flag, t).*vT;
        % Conservatio neqn
        eqs{i} = (op.pv./dt).*(t - t0) + op.Div(vi);
        
        % Composition source terms
        qi = (inj.*comp(:, i) + ~inj.*t(wc)).*qT;
        eqs{i}(wc) = eqs{i}(wc) - qi;
        
        types{i} = 'cell';
    end
    
    problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end

function [vT, qT] = getFlux(model, state)
    vT = sum(state.flux(model.operators.internalConn, :), 2);
    qT = sum(vertcat(state.wellSol.flux), 2);
end