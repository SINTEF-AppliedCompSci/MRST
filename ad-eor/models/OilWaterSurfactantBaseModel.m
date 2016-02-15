classdef OilWaterSurfactantBaseModel < ReservoirModel 

   methods
      
      function model = OilWaterSurfactantBaseModel(G, rock, fluid, varargin)

         model = model@ReservoirModel(G, rock, fluid, varargin{:});
         model = model.setupOperators(G, rock, varargin{:});
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

   end
end
