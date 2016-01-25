classdef OilWaterSurfactantModel < TwoPhaseOilWaterModel
% Oil/water/Surfactant system
% This model is a two phase oil/water model, extended with the surfactant
% component in addition.

   properties
      % Surfactant present
      surfactant

   end

   methods
      function model = OilWaterSurfactantModel(G, rock, fluid, varargin)

         model = model@TwoPhaseOilWaterModel(G, rock, fluid);

         % This is the model parameters for oil/water/surfactant
         model.surfactant = true;

         model.wellVarNames = {'qWs', 'qOs', 'qWSft', 'bhp'};

         model = merge_options(model, varargin{:});

         model = model.setupOperators(G, rock, varargin{:});

      end

      function [problem, state] = getEquations(model, state0, state, ...
                                                      dt, drivingForces, varargin)
         [problem, state] = equationsOilWaterSurfactant(state0, state, ...
                                                        model, dt, drivingForces, varargin{:});
      end

      function [state, report] = updateState(model, state, problem, ...
                                             dx, drivingForces)
         [state, report] = updateState@TwoPhaseOilWaterModel(model, ...
                                                           state, problem,  dx, drivingForces);

         if model.surfactant
            c = model.getProp(state, 'surfactant');
            state = model.setProp(state, 'surfactant', max(c, 0) );
         end
      end

      function model = setupOperators(model, G, rock, varargin)
         model.operators.veloc = computeVelocTPFA(G, model.operators.internalConn);
      end

      function varargout = evaluateRelPerm(model, sat, varargin)
         error('function evaluateRelPerm is not implemented for surfactant model')
      end


      function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
         [state, report] = updateAfterConvergence@TwoPhaseOilWaterModel(model, state0, state, ...
                                                           dt, drivingForces);
         % Nothing extra to be done here for the surfactant case.
      end


      function [fn, index] = getVariableField(model, name)
      % Get the index/name mapping for the model (such as where
      % pressure or water saturation is located in state)
         switch(lower(name))
           case {'surfactant'}
             index = 1;
             fn = 'c';
           otherwise
             [fn, index] = getVariableField@TwoPhaseOilWaterModel(...
                model, name);
         end
      end

      function scaling = getScalingFactorsCPR(model, problem, names)
         nNames = numel(names);

         scaling = cell(nNames, 1);
         handled = false(nNames, 1);

         for iter = 1:nNames
            name = lower(names{iter});
            switch name
              case 'surfactant'
                s = 0;
              otherwise
                continue
            end
            sub = strcmpi(problem.equationNames, name);

            scaling{iter} = s;
            handled(sub) = true;
         end
         if ~all(handled)
            % Get rest of scaling factors
            other = getScalingFactorsCPR@ThreePhaseBlackOilModel(model, problem, names(~handled));
            [scaling{~handled}] = other{:};
         end
      end

   end
end
