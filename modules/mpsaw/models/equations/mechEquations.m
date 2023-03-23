function [eqs, names, types, state] = mechEquations(model, state0, state, dt, drivingForces, varargin)
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
