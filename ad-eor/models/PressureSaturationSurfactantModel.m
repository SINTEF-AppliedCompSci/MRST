classdef PressureSaturationSurfactantModel < OilWaterSurfactantBaseModel

   methods
      function model = PressureSaturationSurfactantModel(G, rock, fluid, varargin)
         model = model@OilWaterSurfactantBaseModel(G, rock, fluid, varargin{:});
      end

      function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
         [problem, state] = equationsOilWaterSurfactant(state0, state, model, dt, drivingForces, ...
                                                        'assembleOnlyOWEquation', true, ...
                                                        varargin{:});
      end


   end
end
