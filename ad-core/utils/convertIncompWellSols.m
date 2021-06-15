function wellSols = convertIncompWellSols(W, states, incompFluid)
%Convert wellSols from incomp module to format used in ad-core/ad-blackoil
%
% SYNOPSIS:
%   wellSols = convertIncompWellSols(W, states, fluid)
%
% DESCRIPTION:
%   The solvers in the `incomp` module uses a different wellSol format than
%   the `ad-core` style wellSols. This function converts a set of states
%   into wellSols suitable for routines that were designed to work with the
%   `ad-core` style wellSols. Specifically, this enables the use of
%   getWellOutput and plotWellSols with solutions from the incompressible
%   solvers not based on AD.
%
% REQUIRED PARAMETERS:
%   W           - Wells used for simulation
%
%   states      - Nstep long struct array of all simulation states.
%
%   incompFluid - Fluid model used to compute the simulation states.
%
% RETURNS:
%   wellSols    - Nstep long cell array of the same format as output by
%                 e.g. `simulateScheduleAD`
%
% SEE ALSO:
%   `plotWellSols`, `getWellOutput`

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
nstep = numel(states);
wellSols = cell(nstep, 1);
phases = {'W', 'O'};

for i = 1:nstep
    ws = [];
    state = states(i);
    mu = incompFluid.properties(state);
    s         = incompFluid.saturation(state);
    kr        = incompFluid.relperm(s, state);
    
    mob    = bsxfun(@rdivide, kr, mu);
    totmob = sum(mob, 2);
    
    for j = 1:numel(W)
        wold = states(i).wellSol(j);
        
        w.name = W(j).name;
        
        influx = min(wold.flux, 0);
        outflux = max(wold.flux, 0);

        wc = W(j).cells;
        f = bsxfun(@rdivide, mob(wc, :), totmob(wc));
        for k = 1:2
            fn = ['q', phases{k}, 's'];
            w.(fn) = sum(influx.*f(:, k)) + sum(outflux*W(j).compi(k));
        end
        w.bhp = wold.pressure;
        ws = [ws; w];
    end
    wellSols{i} = ws;
end
end
