classdef ExplicitSurfactantModel1 < OilWaterSurfactantModel
% Oil/water/Surfactant system
% This model is a two phase oil/water model, extended with the surfactant
% component in addition.

   methods
      function model = ExplicitSurfactantModel1(G, rock, fluid, varargin)

         model = model@OilWaterSurfactantModel(G, rock, fluid);

         model.wellVarNames = {'qWSft', 'bhp'};
         model = merge_options(model, varargin{:});


      end

      function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsOilWaterSurfactant(state0, state, model, dt, drivingForces, ...
                                                           'assembleOnlyExplicitConcentrationEquation', ...
                                                           true, varargin{:});
      end

      function [state, report] = updateState(model, state, problem, dx, drivingForces)

         [state, report] = ReservoirModel.updateStateWithWells(model, state, problem,  dx, ...
                                                           drivingForces);
         % cap the concentration (only if implicit solver for concentration)
         c = model.getProp(state, 'surfactant');
         state = model.setProp(state, 'surfactant', max(c, 0) );

      end

      function varargout = evaluateRelPerm(model, sat, varargin)
         error('function evaluateRelPerm is not implemented for surfactant model')
      end

      function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)

         [state, report] = updateAfterConvergence@TwoPhaseOilWaterModel(model, state0, state, dt, ...
                                                           drivingForces);
         % Note that the adsorption values are not updated.

      end


   end
end
