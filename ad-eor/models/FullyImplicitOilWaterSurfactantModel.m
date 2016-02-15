classdef FullyImplicitOilWaterSurfactantModel < OilWaterSurfactantBaseModel
% Oil/water/Surfactant system
% This model is a two phase oil/water model, extended with the surfactant
% component in addition.


   methods

      function model = FullyImplicitOilWaterSurfactantModel(G, rock, fluid, varargin)
         model = model@OilWaterSurfactantBaseModel(G, rock, fluid, varargin{:});
      end

      function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsOilWaterSurfactant(state0, state, model, dt, drivingForces, ...
                                                           varargin{:});
      end

      function [state, report] = updateState(model, state, problem, dx, drivingForces)
         [state, report] = updateState@OilWaterSurfactantBaseModel(model, state, problem,  dx, ...
                                                           drivingForces);
         % cap the concentration (only if implicit solver for concentration)
         c = model.getProp(state, 'surfactant');
         state = model.setProp(state, 'surfactant', max(c, 0) );

      end

      function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
         [state, report] = updateAfterConvergence@TwoPhaseOilWaterModel(model, state0, state, dt, ...
                                                           drivingForces);
         state = updateAdsorption(state0, state, model);

      end

   end
end
