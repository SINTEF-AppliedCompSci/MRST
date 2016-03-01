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

         model.pressureSaturationSurfactantSolver = NonLinearSolver;
         model.explicitConcentrationSolver = NonLinearSolver;

      end

      
      function [state, report] = stepFunction(model, state, state0, dt, ...
                                             drivingForces, linsolve, nonlinsolve,...
                                              iteration, varargin)

         press_sat_model = model.pressureSaturationSurfactantModel;
         conc_model = model.explicitConcentrationModel;

         press_sat_solver = model.pressureSaturationSurfactantSolver;
         conc_solver = model.explicitConcentrationSolver;
         
         forceArg = getDrivingForces(press_sat_model, drivingForces);
         
         [state, press_sat_report] = press_sat_solver.solveTimestep(state0, dt, press_sat_model, ...
                                                           'initialGuess', state, forceArg{:});

         if press_sat_report.Converged
            % Solve for concentration
            [new_state, conc_report] = conc_solver.solveTimestep(state, dts(subiter), conc_model, ...
                                                              forceArg{:});
            new_state = updateAdsorption(state, new_state, model);
            state = new_state;
         else
            conc_report.Converged = false;
         end
         
         converged = press_sat_report.Converged && conc_report.Converged;
         failure = press_sat_report.failure && conc_report.failure;

         report = model.makeStepReport(...
            'Converged',       converged, ...
            'failure', failure ...
            );
         
         report.PressureSaturationSolver =  press_sat_report;
         report.ConcentrationSolver= conc_report;
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