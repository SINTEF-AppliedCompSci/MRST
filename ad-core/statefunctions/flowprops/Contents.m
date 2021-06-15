% FLOWPROPS
%
% Files
%   BaseRelativePermeability          - Phase relative permeability (per cell)
%   BlackOilCapillaryPressure         - Implementation of black-oil type capillary pressure
%   BlackOilDensity                   - 
%   BlackOilPoreVolume                - Effective pore-volume after pressure-multiplier
%   BlackOilPressureReductionFactors  - Component weighting factors used to form a pressure equation
%   BlackOilShrinkageFactors          - Shrinkage factors for black-oil
%   BlackOilViscosity                 - Black-oil style viscosity functions that account for rs and Rv
%   ComponentMobility                 - Class implementing the mobility for each component and phase
%   ComponentPhaseDensity             - Component density in each cell for each phase
%   ComponentPhaseMass                - Mass of each component, in each phase.
%   ComponentProperty                 - Virtual class for properties that get their values from the model's
%   ComponentTotalMass                - The total mass of each component, given per cell
%   DensityDerivedShrinkageFactors    - Simple "shrinkage factors" which are just phase densities divided by
%   Mobility                          - Mobility for each phase. Normally rel. perm. divided by viscosity.
%   NumericalPressureReductionFactors - 
%   PhasePressures                    - Pressure for each phase. Will always return one value for each phase,
%   PoreVolume                        - Static pore-volume taken from state
%   PressureReductionFactors          - Component weighting factors used to form a pressure equation -
%   RsMax                             - Maximum Rs (dissolved gas-oil ratio)
%   RvMax                             - Maximum Rv (vaporized oil-gas ratio)
%   ShrinkageFactors                  - Shrinkage factors that depend only on pressure
%   Viscosity                         - 

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
