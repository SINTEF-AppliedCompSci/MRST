function state = Advance(model, dt, scheduleRow)      
    state   = model.state;    
    fluid   = model.fluid; 
    psolve  = model.solver.psolve;
    tsolve  = model.solver.tsolve;
    type    = model.experiment.process.type;
    g       = scheduleRow(3);
    
    if(strcmp(type,'CENT'))
        gravity('on','x', -g)
    else
        gravity off
    end
    % solve transport
    state = tsolve(state, dt, fluid);
    % Check for inconsistent saturations
    s = state.s(:,1);
    assert(max(s) < 1+eps && min(s) > -eps);
    % update pressure solution
    state = psolve(state, fluid);    
end