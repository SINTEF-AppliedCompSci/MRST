classdef ImplicitExplicitOilWaterSurfactantModel < OilWaterSurfactantBaseModel
%
%
% SYNOPSIS:
%   model = ImplicitExplicitOilWaterSurfactantModel(G, rock, fluid, varargin)
%
% DESCRIPTION: Implicit-explicit model for a oil water system with
% surfactant. The implicit step consists of solving the oil-water system
% implicitly with a fixed concentration of surfactant. Then, the transport
% equation for the surfactant are solved explicitly, for fixed pressure and
% saturation. A description of the surfactant model that is implemented can be
% found in ad-eor/docs directory.
%
% PARAMETERS:
%   G        - Grid
%   rock     - Rock structure
%   fluid    - Fluid structure
%   varargin - optional parameter
%
% RETURNS:
%   class instance
%
% EXAMPLE:
%
% SEE ALSO: OilWaterSurfactantBaseModel, FullyImplicitOilWaterSurfactantModel,
% pressureSaturationSurfactantModel, explicitConcentrationModel
%

    properties

        % Instance of the model which is use for solving implicitly the pressure and
        % saturation equation, only, with given surfactant concentration
        pressureSaturationSurfactantModel;

        % Instance of the model which is used for solving the concentration equation
        % explicitly.
        explicitConcentrationModel;

        % Solver to be used for the pressure-saturation equations.
        pressureSaturationSurfactantSolver;

        % Solver to be used for the concentration equations.
        explicitConcentrationSolver;

    end

    methods
        function model = ImplicitExplicitOilWaterSurfactantModel(G, rock, fluid, varargin)

            model = model@OilWaterSurfactantBaseModel(G, rock, fluid, varargin{:});
            model = merge_options(model, varargin{:});

            % Setup the properties specific to the model, see comments in
            % properties block.

            model.pressureSaturationSurfactantModel = PressureSaturationSurfactantModel(G, rock, fluid, varargin{:});
            model.explicitConcentrationModel = ExplicitConcentrationModel(G, rock, fluid, varargin{:});

            model.pressureSaturationSurfactantSolver = NonLinearSolver();
            model.explicitConcentrationSolver = NonLinearSolver();

            model = merge_options(model, varargin{:});

        end

        function [state, report] = stepFunction(model, state, state0, dt, ...
                                                drivingForces, linsolve, nonlinsolve,...
                                                iteration, varargin)

            press_sat_model = model.pressureSaturationSurfactantModel;
            conc_model = model.explicitConcentrationModel;

            press_sat_solver = model.pressureSaturationSurfactantSolver;
            conc_solver = model.explicitConcentrationSolver;

            % Solve the pressure - saturation equation for the time step dt.
            [state, press_sat_converged] = solveMinistep(press_sat_solver, press_sat_model, state, state0, dt, drivingForces);

            if press_sat_converged
                % If the pressure-saturation equations have been solved
                % successfully, solve the concentration equation.
                [new_state, conc_converged] = solveMinistep(conc_solver, conc_model, state, state0, ...
                                                            dt, drivingForces);
                if conc_converged
                    new_state = updateAdsorption(state, new_state, model);
                    state = new_state;
                end
            end
            converged = press_sat_converged && conc_converged;
            report = model.makeStepReport('Converged', converged);

        end

    end
end

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
