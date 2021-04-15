% ADJOINT Add-on module for adjoint-based production optimisation
%
% Files
%   addAdjointWellFields            - SYNOPSIS:
%   adjointFluidFields              - Extend fluid functionality with fields needed in (2nd order) adjoint imp.
%   assembleWellSystem              - Generate pressure linear system components for wells.
%   computeAdjointRHS               - Compute adjoint 'pressure' rhs
%   computeGradient                 - compute gradient for control variables and project according to
%   computeNumericalGradient        - compute numerical gradient
%   controls2RHS                    - Create mappings A_N, b_N, A_D, b_D such that
%   controls2Wells                  - Create mappings A_N, b_N, A_D, b_D such that
%   generateUpstreamTransportMatrix - generateUpstreamTransportMatrix for use in saturation solver
%   initControls                    - initControls -- Initialize control structure based on well schedule
%   initSchedule                    - initSchedule -- Initialize schedule structure based on well W.
%   optimizeObjective               - optimizeObjective -- Run whole optimization proccess using ad-hoc line
%   projectGradient                 - Project gradient according to linear input constraints. Handles box-constraints and
%   runAdjoint                      - runAdjoint -- Run adjoint simulation based on simRes and schedule.
%   runSchedule                     - runSchedule -- Run simulation based on schedule.
%   solveAdjointPressureSystem      - Find current time step (search for empty slots in adjRes)
%   solveAdjointTransportSystem     - Find current time step (search for empty slots in adjRes)
%   updateSchedule                  - Update schedule based on controls
%   updateWells                     - Update wells based on schedule time step
%   dispControls                    - Display control values.
%   dispSchedule                    - Display schedule values
%   lineSearchAgr                   - Run agressive line search based on given gradient.
%   solveIncompFlowLocal            - Local version of solveIncompFlow for use with adjoint module. Local

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
