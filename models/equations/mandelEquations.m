function [eqs, names, types, state] = mandelEquations(model, state0, state, dt, drivingForces, varargin)
            
% The mechanical part of the assembly takes the form
% 
% | A   -D   0  0 |   |u  |   |  0 |
% | D'   0   0  R | * |lm | = |  0 |
% | 0    0   R' 0 |   |vd |   | -F |
%
% where u : displacement
%      lm : Lagrangian multiplier for displacement constraint at boundary
%      vd : Vertical displacement (scalar value) of the top plate
%      F  : Total force exterted at the top
%
    
    G     = model.G;
    fluid = model.fluid;
    op    = model.operators;
    
    pv         = op.pv;
    divKgradop = op.divKgradop;
    divuop     = op.divuop;
    momentop   = op.momentop;
    mechdirop  = op.mechDirichletop;
    fluiddirop = op.fluidDirichletop;
    R = op.R;
    
    u  = model.getProp(state, 'u');
    p  = model.getProp(state, 'p');
    lm = model.getProp(state, 'lambdamech');
    lf = model.getProp(state, 'lambdafluid');
    vd = model.getProp(state, 'vertdisp');
    avgtopforce = model.getProp(state, 'avgtopforce');
    
    u0  = model.getProp(state0, 'u');
    p0  = model.getProp(state0, 'p');
    lm0 = model.getProp(state0, 'lambdamech');
    
    c = fluid.c;
    fac =  1 / (1e6 * mean(G.cells.volumes));
    
    divu  = model.getProp(state, 'Dilatation');
    divu0 = model.getProp(state0, 'Dilatation');
    
    eqs{1} = fac*momentop(u, p, lm);
    eqs{2} = 1/dt.*(pv.*c.*(p - p0) + (divu - divu0)) + divKgradop(p, lf);
    eqs{3} = fac*(mechdirop(u, p, lm) - R*vd);
    eqs{4} = fluiddirop(p, lf);
    eqs{5} = R'*lm + avgtopforce;
    
    names = {'momentum', 'mass', 'bcmech', 'bcfluid', 'vertdisp'};
    types = {'cellcol', 'cell', 'bcmech', 'bcfluid', 'vertdisp'};
    
end

