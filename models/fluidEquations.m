function [eqs, names, types, state] = fluidEquations(model, state0, state, dt, drivingForces, varargin)
            
    G   = model.G;
    op  = model.operators;
    B   = op.B;
    rhs = op.rhs;
    
    c = model.fluid.c;
    
    p      = model.getProp(state, 'p');
    p0     = model.getProp(state0, 'p');
    lambda = model.getProp(state, 'lambdafluid');
    pv     = model.getProp(state, 'PoreVolume');

    fac = max(rhs{2});
    
    accterm = c*pv.*(p - p0);
    eqs{1} = 1/dt*accterm + (B{1,1}*p + B{1,2}*lambda - rhs{1});
    eqs{2} = 1/fac*(B{2,1}*p + B{2,2}*lambda - rhs{2});
    
    eqs = vertcat(eqs{:});
            
    names = {'fluid'};
    types = {'mixed'};
end

