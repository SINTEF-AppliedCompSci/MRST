classdef CO2Model < TwoPhaseOilWaterModel
% Function to implement the two-phase flow model to asses the co2 leakage.
% 
% This function is modified from a file in The MATLAB Reservoir Simulation
% Toolbox (MRST), see
%   mrst/modules/ad-eor/models/OilWaterPolymerModel.m
%
% We refer to that function for a complete commented version of the file. 

%{ 
Partial copyright 2009-2021, SINTEF Digital, Mathematics & Cybernetics.
Partial copyright 2021, NORCE Norwegian Research Centre AS, Computational 
Geosciences and Modeling.

This file is part of the ad-micp module.

ad-micp is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ad-micp is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.  

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%}

    methods
        function model = CO2Model(G, rock, fluid, varargin)
            model = model@TwoPhaseOilWaterModel(G, rock, fluid);
            model = merge_options(model, varargin{:});
        end

        function [problem, state] = getEquations(model, state0, state, ...
                                               dt, drivingForces, varargin)
            [problem, state] = equationsCO2(state0, state, ...
                                    model, dt, drivingForces, varargin{:});
        end

        function state = validateState(model, state)
            state = validateState@TwoPhaseOilWaterModel(model, state);
        end

        function [state, report] = updateState(model, state, problem, ...
                                                         dx, drivingForces)
            [state, report] = updateState@TwoPhaseOilWaterModel(model, ...
                                       state, problem,  dx, drivingForces);
        end

        function [state, report] = updateAfterConvergence(model, ...
                                          state0, state, dt, drivingForces)
            [state, report] = ...
                    updateAfterConvergence@TwoPhaseOilWaterModel(model, ...
                                         state0, state, dt, drivingForces);
        end

        function [fn, index] = getVariableField(model, name, varargin)
                    [fn, index] = ...
                          getVariableField@TwoPhaseOilWaterModel(model, ...
                                                        name, varargin{:});
        end

        function names = getComponentNames(model)
            names = getComponentNames@TwoPhaseOilWaterModel(model);
        end

        function [compEqs, compSrc, eqNames, wellSol] = ...
                       getExtraWellContributions(model, well, wellSol0, ...
                      wellSol, q_s, bh, packed, qMass, qVol, dt, iteration)
            [compEqs, compSrc, eqNames, wellSol] = ...
                 getExtraWellContributions@TwoPhaseOilWaterModel(model, ...
                 well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, ...
                                                            dt, iteration);
        end
    end
end