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
         press_sat_ok = press_sat_report.Converged;
         
         if press_sat_ok
            % Solve for concentration
            dts = splitTime(dt, 0.1*day);
            dts = dt;
            stateM = state0;
            conc_ok = true;
            subiter = 1;
            while conc_ok & (subiter <= numel(dts))
               [state, conc_report] = conc_solver.solveTimestep(stateM, dts(subiter), conc_model, ...
                                                                'initialGuess', state, forceArg{:});
               state = updateAdsorption(stateM, state, model);
               stateM = state;
               conc_ok = conc_report.Converged;
               subiter = subiter + 1;
            end
            state.SURFADS = double(state.ads);
         else
            conc_ok = false;
            conc_report = [];
         end

         
         converged = press_sat_ok && conc_ok;

         if ~press_sat_ok
            FailureMsg = 'Pressure failed to converge!';
         else
            FailureMsg = '';
         end

         report = model.makeStepReport(...
            'Failure',        ~press_sat_ok, ...
            'Converged',       converged, ...
            'FailureMsg',      FailureMsg ...
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