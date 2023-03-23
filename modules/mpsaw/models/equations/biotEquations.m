function [eqs, names, types, state] = biotEquations(model, state0, state, dt, drivingForces, varargin)
%Undocumented Utility Function

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
