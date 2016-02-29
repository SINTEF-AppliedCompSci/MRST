classdef PressureSaturationSurfactantModel < OilWaterSurfactantBaseModel

   methods
      function model = PressureSaturationSurfactantModel(G, rock, fluid, varargin)
         model = model@OilWaterSurfactantBaseModel(G, rock, fluid, varargin{:});
      end

      function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
         [problem, state] = equationsPressureSaturationForOilWaterSurfactant(state0, state, model, ...
                                                           dt, drivingForces, varargin{:});
      end


   end
end
