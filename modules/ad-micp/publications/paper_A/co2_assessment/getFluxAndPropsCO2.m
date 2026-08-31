function [co2Flux, co2Mobility] = ...
        getFluxAndPropsCO2( ...
        model, co2Pressure, co2RelativePermeability, ...
        transmissibility, gravityGradient, varargin)
% Compute CO2 flux and cell mobility.
%
% The non-wetting CO2 phase is represented internally by the oil phase of
% MRST`s TwoPhaseOilWaterModel.
%
% This function is modified from a file in the MATLAB Reservoir Simulation
% Toolbox (MRST), see
%   mrst/modules/ad-eor/utils/getFluxAndPropsWaterPolymer_BO.m
%
% We refer to that function for a complete commented version of the file.

%{
Partial copyright 2009-2021, SINTEF Digital, Mathematics & Cybernetics.
Partial copyright 2021-2026, NORCE Research AS, Computational Geosciences and
Modeling.

This file is part of the ad-micp module.

ad-micp is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ad-micp is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file. If not, see <http://www.gnu.org/licenses/>.
%}
    fluid = model.fluid;
    operators = model.operators;

    co2PotentialDifference = ...
        operators.Grad(co2Pressure) - ...
        fluid.rhoOS .* gravityGradient;

    upstreamCO2 = ...
        value(co2PotentialDifference) <= 0;

    [faceRelativePermeability, cellRelativePermeability] = ...
        operators.splitFaceCellValue( ...
        operators, ...
        upstreamCO2, ...
        co2RelativePermeability);

    co2Mobility = ...
        cellRelativePermeability ./ fluid.muO;

    faceMobility = ...
        faceRelativePermeability ./ fluid.muO;

    co2Flux = ...
        -faceMobility .* ...
        transmissibility .* ...
        co2PotentialDifference;
end
