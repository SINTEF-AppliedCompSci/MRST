% UTILS
%
% Files
%   calculateHydrocarbonsFromStatusBO - Compute solution variables for the gas/oil/rs/rv-variable in black-oil
%   computeFlashBlackOil              - Compute flash for a black-oil model with disgas/vapoil
%   computeInitAquifer                - Undocumented Utility Function
%   equationsBlackOil                 - Generate linearized problem for the black-oil equations
%   equationsOilWater                 - Generate linearized problem for the two-phase oil-water model
%   equationsWater                    - Generate linearized problem for the single-phase water model
%   equationsWaterThermal             - Get linearized problem for single phase water system with black oil-style
%   getbG_BO                          - Utility function for evaluating the reciprocal gas formation volume
%   getbO_BO                          - Utility function for evaluating the reciprocal oil formation volume
%   getCapillaryPressureBO            - Compute oil-water and oil-gas capillary pressures
%   getCellStatusVO                   - Get status flags for each cell in a black-oil model
%   getDeckEGG                        - Get the parsed deck for the EGG model
%   getFluxAndPropsGas_BO             - Get flux and properties for the gas phase for a black-oil problem
%   getFluxAndPropsOil_BO             - Get flux and properties for the oil phase for a black-oil problem
%   getFluxAndPropsWater_BO           - Get flux and properties for the water phase for a black-oil problem
%   getPolymerShearMultiplier         - Compute the flux multiplier due to polymer shear thinning/thickening
%   getWellPolymer                    - Undocumented Utility Function
%   updateStateBlackOilGeneric        - Generic update function for blackoil-like models

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
