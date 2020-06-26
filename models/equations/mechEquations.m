function [eqs, names, types, state] = mechEquations(model, state0, state, dt, drivingForces, varargin)
            
    G   = model.G;
    op  = model.operators;
    B   = op.B;
    rhs = op.rhs;
    
    u         = model.getProp(state, 'u');
    lambda    = model.getProp(state, 'lambdamech');
    
    eqs{1} = B{1,1}*u + B{1,2}*lambda - rhs{1};
    eqs{2} = B{2,1}*u + B{2,2}*lambda - rhs{2};
            
    names = {'momentum', 'bc'};
    types = {'cells', 'bc'};
    
end

