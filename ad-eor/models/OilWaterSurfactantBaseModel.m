classdef OilWaterSurfactantBaseModel < TwoPhaseOilWaterModel

   properties
      surfactant
   end

   methods

      function model = OilWaterSurfactantBaseModel(G, rock, fluid, varargin)

         model = model@TwoPhaseOilWaterModel(G, rock, fluid, varargin{:});
         model = model.setupOperators(G, rock, varargin{:});
         model.surfactant = true;
         model.wellVarNames = {'qWs', 'qOs', 'qWSft', 'bhp'};
         model = merge_options(model, varargin{:});

      end

      function model = setupOperators(model, G, rock, varargin)
         model.operators.veloc = computeVelocTPFA(G, model.operators.internalConn);
         model.operators.sqVeloc = computeSqVelocTPFA(G, model.operators.internalConn);
      end


      function varargout = evaluateRelPerm(model, sat, varargin)
         error('function evaluateRelPerm is not implemented for surfactant model')
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

      function state = storeSurfData(model, state, s, c, Nc, sigma)
         state.SWAT    = double(s);
         state.SURFACT = double(c);
         state.SURFCNM = log(double(Nc))/log(10);
         state.SURFST  = double(sigma);
         % state.SURFADS = double(ads);
      end

   end
end
