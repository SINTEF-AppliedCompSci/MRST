function [eqs, names, types, state] = biotEquations(model, state0, state, dt, drivingForces, varargin)
            
    G     = model.G;
    fluid = model.fluid;
    op    = model.operators;

    pv         = op.pv;
    divKgradop = op.divKgradop;
    divuop     = op.divuop;
    momentop   = op.momentop;
    mechdirop  = op.mechDirichletop;
    fluiddirop = op.fluidDirichletop;
    
    
    u  = model.getProp(state, 'u');
    p  = model.getProp(state, 'p');
    p0 = model.getProp(state0, 'p');
    lm = model.getProp(state, 'lambdamech');
    lf = model.getProp(state, 'lambdafluid');
    
    c = fluid.c;
    
    eqs{1} = momentop(u, p, lm);
    eqs{2} = 1/dt*(c*pv.*(p - p0) + divuop(u, p, lm)) + divKgradop(p, lf);
    eqs{3} = mechdirop(u, p, lm);
    eqs{4} = fluiddirop(p, lf);
    
    names = {'momentum', 'mass', 'bcmech', 'bcfluid'};
    types = {'cellcol', 'cell', 'bcmech', 'bcfluid'};
    
end

