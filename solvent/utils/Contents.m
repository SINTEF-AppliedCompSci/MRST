% UTILS
%
% Files
%   addSolventProperties           - Add solvent pseudo-component and properties to fluid.
%   computeFlashBlackOilSolvent    - Compute flash for a black-oil solvent model with disgas/vapoil
%   computeRelPermSolvent          - Calulates effective relative permeabilities.
%   computeResidualSaturations     - Calculate effective residual saturations
%   computeViscositiesAndDensities - Calculates effective viscosities and densities using Todd-Longstaff
%   getFluxAndPropsSolvent         - Flux and other properties for the black-oil solvent model equations.
%   getPlotAfterStepSolvent        - Get a function that allows for dynamic plotting in `simulateScheduleAD`
%   makeWAGschedule                - Make schedule for water-alternating gas (WAG) injection.
%   plotSolventFluidProps          - Terneary plots of HC components in black-oil solvent model.

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
