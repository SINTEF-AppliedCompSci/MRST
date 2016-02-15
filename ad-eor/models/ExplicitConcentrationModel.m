classdef ExplicitConcentrationModel < OilWaterSurfactantBaseModel

   methods
      function model = ExplicitConcentrationModel(G, rock, fluid, varargin)
         model = model@OilWaterSurfactantBaseModel(G, rock, fluid, varargin{:});
         model.wellVarNames = {'qWSft', 'bhp'};
      end


      function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsOilWaterSurfactant(state0, state, model, dt, drivingForces, ...
                                                           'assembleOnlyExplicitConcentrationEquation', ...
                                                           true, varargin{:});
      end

      function [state, report] = updateState(model, state, problem, dx, drivingForces)

         [state, report] = updateState@OilWaterSurfactantBaseModel(model, state, problem,  dx, ...
                                                           drivingForces);
         % cap the concentration (only if implicit solver for concentration)
         c = model.getProp(state, 'surfactant');
         state = model.setProp(state, 'surfactant', max(c, 0) );
            
      end
      
      function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
         [state, report] = updateAfterConvergence@OilWaterSurfactantBaseModel(model, state0, state, ...
                                                           dt, drivingForces);
         state = updateAdsorption(state0, state, model);
      end
      
   end
end
