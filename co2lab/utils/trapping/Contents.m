% TRAPPING
%
% Files
%   downstreamTraps           - Return index of all traps downstream of 'trap_ix'.  Also include
%   findOptimalInjectionPoint - Find the optimal point to inject CO2
%   flattenTraps              - Given output res from trapAnalysis, produce a grid with flat top surfaces in trap areas.
%   maximizeTrapping          - Find the N best injection trees, with optional compensation for overlap
%   trapAnalysis              - Compute and summarize the relevant trap analysis information
%   volumesOfTraps            - Compute volumes of (a subset of) precomputed structural traps

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
