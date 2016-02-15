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
         model.wellVarNames = {'qWs', 'qOs', 'qWSft', 'bhp'};
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
            % Solve concentration
            [state, conc_report] = conc_solver.solveTimestep(state0, dt, conc_model, 'initialGuess', ...
                                                             state, forceArg{:});
            conc_ok = conc_report.Converged;
         else
            conc_ok = false;
            conc_report = [];
         end
         values = press_sat_report.StepReports{end}.NonlinearReport{end}.Residuals;
         
         converged = press_sat_ok && con_ok;
         if converged && ~model.stepFunctionIsLinear
            problem = press_sat_model.getEquations(state0, state, dt, drivingForces, 'resOnly', true, 'iteration', inf);
            % Is the pressure still converged when accounting for the
            % mobilities?
            [~, values] = press_sat_model.checkConvergence(problem);
            if model.outerCheckWellConvergence
               converged = all(values < model.outerTolerance);
            else
               converged = values(1) < model.outerTolerance;
            end
            converged = converged || iteration > model.maxOuterIterations;
         end
         if ~press_sat_ok
            FailureMsg = 'Pressure failed to converge!';
         else
            FailureMsg = '';
         end

         
         report = model.makeStepReport(...
            'Failure',        ~pressure_ok, ...
            'Converged',       converged, ...
            'FailureMsg',      FailureMsg, ...
            'Residuals',       values ...
            );
         
         report.press_sat_solver =  press_sat_report;
         report.conc_solver= conc_solver;
      end


      function state = storeSurfData(model, state, s, c, Nc, sigma)
         state.SWAT    = double(s);
         state.SURFACT = double(c);
         state.SURFCNM = log(double(Nc))/log(10);
         state.SURFST  = double(sigma);
         % state.SURFADS = double(ads);
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