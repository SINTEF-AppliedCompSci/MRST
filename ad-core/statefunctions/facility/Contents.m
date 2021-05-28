% FACILITY
%
% Files
%   FacilityWellMapping                      - Generate a struct containing a bunch of useful mappings for a given
%   InjectionSurfaceDensity                  - Get injection surface density
%   PerforationComponentPhaseDensity         - Component density to used for each well connection
%   PerforationMobility                      - Mobility in each perforated cell of a well
%   PerforationPressureGradient              - Calculate the pressure difference between the reservoir cells and the
%   WellComponentPhaseFlux                   - Component flux in each phase for wells
%   WellComponentTotalFlux                   - Component total flux for wells (with treatment for cross-flow)
%   WellComponentTotalFluxDensityMix         - Component total flux for wells (with treatment for cross-flow)
%   WellComponentTotalVolumeBalanceCrossflow - Component total flux for wells (with treatment for cross-flow)
%   WellIndex                                - Get the well index (productivity or injectivity index, equivialent to
%   WellPhaseFlux                            - Get phase-flux between well-bore and reservoir

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
