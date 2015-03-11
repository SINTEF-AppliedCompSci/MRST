classdef CO2VEBlackOilTypeModel < ReservoirModel
   
   % ============================= Class properties ==========================
   properties
      disgas;  % if 'true', include dissolution effect
      disrate; % Dissolution rate (if dissolution is active).  A value of 0
               % means instantaneous dissolution.
   end
      
   % ============================== Public methods ===========================
   methods

      %% Constructor
      function model = CO2VEBlackOilTypeModel(Gt, rock2D, fluid, varargin)
         
         opt.disgas = false; % default: do not include dissolution effect
         opt.disrate = 0;    % default: instantaneous dissolution (if included)
         opt = merge_options(opt, varargin{:});
         
         model@reservoirModel(Gt, rock2d, fluid);
         model.water   = true;
         model.gas     = true;
         model.oil     = false;
         model.disgas  = opt.disgas;
         model.disrate = opt.disrate;
         
      end
      
   % =========================== Private methods ============================
      function model = setupOperators(model, Gt, rock, varargin)
         
         % @@ remove need for 'useNewStandard' parameter passing
         model.operators = ...
             setupSimCompVe(Gt, rock, 'useNewStandard', true, varargin{:}); 
         
         % @@ In case of a height-formulation, operators.pv should be divided
         % by height, i.e. model.operators.pv = model.operators.pv ./ Gt.cells.H
      end
   % ------------------------------------------------------------------------
      function [problem, state] = ...
             getEquations(model, state0, state, dt, drivingForces, varargin)
      
      end
   
   % ------------------------------------------------------------------------
      function [fn, index] = getVariableField(model, name)
   
      end
   % ------------------------------------------------------------------------
   
      function [state, report] = updateState(model, state, problem, dx, drivingForces)
         
      end
end   
