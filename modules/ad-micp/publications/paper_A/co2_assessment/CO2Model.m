classdef CO2Model < TwoPhaseOilWaterModel
% Model for two-phase water and CO2 flow.
%
% This class uses the standard MRST two-phase oil-water infrastructure,
% with the non-wetting oil phase representing CO2.
%
% This class is modified from a file in the MATLAB Reservoir Simulation
% Toolbox (MRST), see
%   mrst/modules/ad-eor/models/OilWaterPolymerModel.m
%
% We refer to that class for a complete commented version of the model
% interface.

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
    methods
        function model = CO2Model(G, rock, fluid, varargin)
            model = model@TwoPhaseOilWaterModel(G, rock, fluid);
            model = merge_options(model, varargin{:});
        end

        function [problem, state] = getEquations( ...
                model, state0, state, dt, drivingForces, varargin)

            [problem, state] = equationsCO2( ...
                state0, ...
                state, ...
                model, ...
                dt, ...
                drivingForces, ...
                varargin{:});
        end
    end
end
