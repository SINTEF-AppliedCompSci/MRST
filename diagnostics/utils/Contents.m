% Files
%   computeFandPhi.m        - Compute flow-capacity/storage-capacity diagram (F,Phi)
%   computeLorenz.m         - Compute the Lorenz coefficient
%   computeSweep.m          - Compute sweep efficiency versus dimensionless time (PVI)
%   computeTOFandTracer.m   - Compute time-of-flight and tracer distribution using finite-volume scheme.
%   computeWellPairs.m      - Compute volumes and fluxes associated with each flux pair
%   expandWellCompletions.m - Pseudo-wells for computation of flow diagnostics for completions
%   expandCoarseWellCompletions.m  - Pseudo-wells for flow diagnostics in coarsened models
%   interactiveDiagnostics.m       - Launch an interactive diagnostics session
%   plotTracerBlend.m              - Plot tracer partition, gray color for multiple tracers present
%   plotWellAllocationComparison.m - Panel comparing well allocation factors for two different models
%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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
