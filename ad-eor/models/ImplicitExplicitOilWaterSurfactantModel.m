classdef ImplicitExplicitOilWaterSurfactantModel < OilWaterSurfactantBaseModel
% Oil/water/Surfactant system
% This model is a two phase oil/water model, extended with the surfactant
% component in addition.

    properties

        pressureSaturationSurfactantModel;
        explicitConcentrationModel;

        pressureSaturationSurfactantSolver;
        explicitConcentrationSolver;

    end

    methods
        function model = ImplicitExplicitOilWaterSurfactantModel(G, rock, fluid, varargin)

            model = model@OilWaterSurfactantBaseModel(G, rock, fluid, varargin{:});
            model = merge_options(model, varargin{:});

            model.pressureSaturationSurfactantModel = PressureSaturationSurfactantModel(G, rock, fluid, varargin{:});
            model.explicitConcentrationModel = ExplicitConcentrationModel(G, rock, fluid, varargin{:});

            model.pressureSaturationSurfactantSolver = NonLinearSolver();
            model.explicitConcentrationSolver = NonLinearSolver();

        end

        function [state, report] = stepFunction(model, state, state0, dt, ...
                                                drivingForces, linsolve, nonlinsolve,...
                                                iteration, varargin)

            press_sat_model = model.pressureSaturationSurfactantModel;
            conc_model = model.explicitConcentrationModel;

            press_sat_solver = model.pressureSaturationSurfactantSolver;
            conc_solver = model.explicitConcentrationSolver;

            [state, press_sat_converged] = solveMinistep(press_sat_solver, press_sat_model, state, state0, dt, drivingForces);

            if press_sat_converged
                % Solve for concentration
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


function dts = splitTime(dt, mini_dt)
    if mini_dt >= dt
        dts = dt;
        return
    else
        n = floor(dt/mini_dt);
        dts = mini_dt*ones(n, 1);
        dts = [dts; dt - n*mini_dt];
        dts = dts(dts>0);
    end
end