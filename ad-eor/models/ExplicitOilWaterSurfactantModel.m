classdef ExplicitOilWaterSurfactantModel < TwoPhaseOilWaterModel
% Oil/water/Surfactant system
% This model is a two phase oil/water model, extended with the surfactant
% component in addition.

   properties
      % Surfactant present
      surfactant;
      explicitConcentration
      % explConcModel
      explConcModel;
   end

   methods
      function model = ExplicitOilWaterSurfactantModel(G, rock, fluid, varargin)

         model = model@TwoPhaseOilWaterModel(G, rock, fluid);

         % This is the model parameters for oil/water/surfactant
         model.surfactant = true;

         model.wellVarNames = {'qWs', 'qOs', 'qWSft', 'bhp'};

         model = merge_options(model, varargin{:});

         model = model.setupOperators(G, rock, varargin{:});
         
         explicitConcentration = true;
         model.explConcModel = ExplicitSurfactantModel(G, rock, fluid, varargin{:});

      end

      function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsOilWaterSurfactant(state0, state, model, dt, drivingForces, 'assembleOnlyOWEquation', true, varargin{:});
      end

      function [state, report] = updateState(model, state, problem, dx, drivingForces)
         [state, report] = updateState@TwoPhaseOilWaterModel(model, state, problem,  dx, ...
                                                           drivingForces);
         if ~model.explicitConcentration
            % cap the concentration (only if implicit solver for concentration)
            c = model.getProp(state, 'surfactant');
            state = model.setProp(state, 'surfactant', max(c, 0) );
         end

      end

      function model = setupOperators(model, G, rock, varargin)
         model.operators.veloc = computeVelocTPFA(G, model.operators.internalConn);
         model.operators.sqVeloc = computeSqVelocTPFA(G, model.operators.internalConn);
      end

      function varargout = evaluateRelPerm(model, sat, varargin)
         error('function evaluateRelPerm is not implemented for surfactant model')
      end


      function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
         [state, report] = updateAfterConvergence@TwoPhaseOilWaterModel(model, state0, state, dt, ...
                                                           drivingForces);

         % Solve the explicit concentration equation.
         solver = NonLinearSolver;
         nIter = 10;
         ddt = dt/nIter;
         stateM = state0;
         for i = 1 : nIter
            [state, converged, failure, its, reports] = solver.solveMinistep(solver, model.explConcModel, ...
                                                              state, stateM, ddt, drivingForces);
            stateM = state;
         end
         state.SURFADS = double(state.ads);
      
      end


      function [fn, index] = getVariableField(model, name)
      % Get the index/name mapping for the model (such as where
      % pressure or water saturation is located in state)
         switch(lower(name))
           case {'ads'} % needed when model.explicitAdsorption
             index = 1;
             fn = 'ads';
           case {'adsmax'} % needed when model.explicitAdsorption
             index = 1;
             fn = 'adsmax';
           case {'surfactant'}
             index = 1;
             fn = 'c';
           case {'surfactantmax'}
             index = 1;
             fn = 'cmax';
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

      function state = storeSurfData(model, state, s, c, Nc, sigma)
         state.SWAT    = double(s);
         state.SURFACT = double(c);
         state.SURFCNM = log(double(Nc))/log(10);
         state.SURFST  = double(sigma);
         % state.SURFADS = double(ads);
      end

   end
end
