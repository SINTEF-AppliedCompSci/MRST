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
    lm = model.getProp(state, 'lambdamech');
    lf = model.getProp(state, 'lambdafluid');
    ef = model.getProp(state, 'extforce');
    
    u0  = model.getProp(state0, 'u');
    p0  = model.getProp(state0, 'p');
    lm0 = model.getProp(state0, 'lambdamech');
    
    c = fluid.c;
    fac =  1 / (1e6 * mean(G.cells.volumes));
    
    divu  = model.getProp(state, 'Dilatation');
    divu0 = model.getProp(state0, 'Dilatation');
    
    % momentum conservation
    eqs{1} = fac*momentop(u, p, lm, ef);
    % mass conservation
    eqs{2} = 1/dt.*(pv.*c.*(p - p0) + (divu - divu0)) + divKgradop(p, lf);
    % bc constraint for mechanics
    eqs{3} = fac*mechdirop(u, p, lm);
    % bc constraint for flow
    eqs{4} = fluiddirop(p, lf);
    
    names = {'momentum', 'mass', 'bcmech', 'bcfluid'};
    types = {'cellcol', 'cell', 'bcmech', 'bcfluid'};
    
end

