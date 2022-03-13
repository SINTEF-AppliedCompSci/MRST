% UTILS
%
% Files
%   averageRock                   - Average version of rock for use in vertical averaging
%   initResSolVE                  - Wrapper for initResSol which adds any extra properties needed by the
%   initResSolVE_s                - Wrapper for initResSol which adds any extra properties needed by the
%   makeReports                   - This function does intermediate processing of simulation data in order to
%   massAtInfinity                - Forecast amount of co2 (in kg) to remain in formation by time infinity.
%   massTrappingDistributionVEADI - Compute the trapping status distribution of CO2 in each cell of a top-surface grid
%   migrateInjection              - Run a simple injection scenario and visualize each time step
%   phaseMassesVEADI              - Compute column masses of undissolved gas, fluid, and dissolved gas.
%   volumesVE                     - SYNOPSIS:

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
