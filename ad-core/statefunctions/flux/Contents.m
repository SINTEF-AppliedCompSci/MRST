% FLUX
%
% Files
%   ComponentPhaseFlux                  - Flux of each component, in each phase
%   ComponentTotalFlux                  - Total flux of all components, summed up over all phases
%   FaceComponentMobility               - Mobility of the component mass for all components in all phases, with
%   FaceMobility                        - Phase mobility, with phase upwind (or alternatively some other
%   GravityPotentialDifference          - Difference in phase potential over a face due to gravity
%   PermeabilityGradientDiscretization  - Class defining the discretization of the permeability
%   PermeabilityPotentialGradient       - The difference in potential over a face, multiplied with a
%   PhaseFlux                           - Phase flux for each phase. The volumetric flux associated with each
%   PhasePotentialDifference            - The potential difference over each interface, for each phase. This is
%   PhasePotentialDifferenceThresholded - 
%   PhaseUpwindFlag                     - Phase-upstream flag for each phase
%   PressureGradient                    - Gradient of phase pressures for internal faces
%   Transmissibility                    - Transmissibility for internal faces. May include an optional
%   TwoPointFluxApproximation           - 

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
