% UTILS
%
% Files
%   assignWellValuesFromControl       - Assign wellSol values when values are set as controls
%   calculateHydrocarbonsFromStatusBO - Compute solution variables for the gas/oil/rs/rv-variable in black-oil
%   computeFlashBlackOil              - Compute flash for a black-oil model with disgas/vapoil
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
%   updateStateBlackOilGeneric        - Generic update function for blackoil-like models

%{
#COPYRIGHT#
%}
