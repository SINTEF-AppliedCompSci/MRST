function solver = SetupSolver(model, dT,scheduleRow)  
    G       = model.grid.G;
    rock    = model.rock;
    T       = model.grid.T;
    state   = model.state;
    fluid   = model.fluid;
    bc      = model.bc;
    verbose = model.verbose;
    type    = model.experiment.process.type;
    g       = scheduleRow(3);
    
    if(strcmp(type,'CENT'))
        gravity('on','x', -g)
    else
        gravity off
    end
    psolve  = @(state, fluid) ...
        incompTPFA(state, G, T, fluid, 'bc', bc);
    tsolve  = @(state, dT, fluid) ...
        implicitTransport(state, G, dT, rock, fluid, ...
        'verbose', verbose, 'bc', bc);
    solver.psolve = psolve;
    solver.tsolve = tsolve;
end