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

%{
Copyright 2020 University of Bergen and SINTEF Digital, Mathematics & Cybernetics.

This file is part of the MPSA-W module for the MATLAB Reservoir Simulation Toolbox (MRST).

The MPSA-W module is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The MPSA-W module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the MPSA-W module.  If not, see <http://www.gnu.org/licenses/>.
%}


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
