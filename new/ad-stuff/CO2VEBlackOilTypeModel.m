classdef CO2VEBlackOilTypeModel < ReservoirModel
   
% ============================= Class properties ==========================
properties
   
   % Equation is chosen based on whether fluid object includes dissolution
   % effects or not 
   equation
   
end
      
% ============================== Public methods ===========================
methods
   
   %% Constructor
   function model = CO2VEBlackOilTypeModel(Gt, rock2D, fluid, varargin)
   
      opt = struct();
      [opt, unparsed] = merge_options(opt, varargin{:});
   
      model@ReservoirModel(Gt, varargin{:});
      model.fluid   = fluid;
      model.water   = true;
      model.gas     = true;
      model.oil     = false;
      model.saturationVarNames = {'sw', 'sg'}; % @@ Design: ideally, we
                                               % should not have to set
                                               % both this _and_ the
                                               % 'water', 'gas' and 'oil'
                                               % flags above. Check with
                                               % maintainer of parent class.
      model.wellVarNames = {'qWs', 'qGs', 'bhp'};
      model.gravity = [0 norm(gravity)];
      
      if isfield(fluid, 'dis_rate')
         % use model equations with dissolution
         model.equation = @equationsWGVEdisgas;
      else
         % use basic model equations (no dissolution)
         model.equation = @equationsWGVEbasic;
      end
      
      model = model.setupOperators(Gt, rock2D);
      
      % This object and its equation does not support temporal variation
      % in temperatures, so if dependence is detected, throw an error
      if nargin(@fluid.bG) > 1
         error(['This model requires that fluid properties are function of ' ...
                'pressure only.']);
      end
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
      
      [problem, state] = model.equation(model         , ...
                                        state0        , ...
                                        state         , ...
                                        dt            , ...
                                        drivingForces , ...
                                        varargin{:});
   end
   
% ------------------------------------------------------------------------

   function [fn, index] = getVariableField(model, name)
      
      switch(lower(name))
        case {'sgmax'}
          index = 1;
          fn = 'sGmax';
        case {'rs'}
          index = 1;
          fn = 'rs';
        otherwise
          [fn, index] = getVariableField@ReservoirModel(model, name);
      end
   end

% ----------------------------------------------------------------------------
function [state, report] = updateState(model, state, problem, dx, drivingForces)

   [state, report] = updateState@ReservoirModel(model, state, problem, dx, ...
                                                drivingForces);
   
   if isfield(model.fluid, 'dis_rate')
      % The model includes dissolution
      if model.fluid.dis_rate > 0
         % rate-driven dissolution
         
         f           = model.fluid;
         sg          = state.s(:,2);
         sg          = min(1, max(0,sg));
         state.s     = [1-sg, sg];    
         state.sGmax = min(1,state.sGmax);
         state.sGmax = max(0,state.sGmax);
         state.sGmax = max(state.sGmax,sg);
         min_rs      = minRs(state.pressure,state.s(:,2),state.sGmax,f,model.G);
         min_rs      = min_rs./state.s(:,1);
         state.rs    = max(min_rs,state.rs);
         state.rs    = min(state.rs,f.rsSat(state.pressure));         
      else
         % instantaneous dissolution
         diff = 1e-3; % @@ magic constant - necessary for convergence in some cases
         state.rs = min(state.rs, model.fluid.rsSat(state.pressure) + diff);
      end
   end
end


   
% ------------------------------------------------------------------------
   
   function [state, report] = ...
       updateAfterConvergence(model, state0, state, dt, drivingForces)
      
      % Here, we update the hysteresis variable 'sGmax'.  If the residual
      % saturation of gas is 0 (i.e. model.fluid.residuals(2) == 0),
      % keeping track of 'sGmax' is not strictly necessary for
      % computation, but it may still be useful information for
      % interpretation, and to simplify program logic we compute it at all
      % times.
      report = []; % not used
      if isfield(model.fluid, 'dis_rate')
         return; % sGmax already updated along with other ADI variables
      end
         
      sGmax0 = model.getProp(state0, 'sGmax');
      sG     = model.getProp(state, 'sg');
      
      state = model.setProp(state, 'sGmax', max(sG, sGmax0));
   end

% ----------------------------------------------------------------------------

   function gdz = getGravityGradient(model)
      s  = model.operators;
      Gt = model.G;
      g  = norm(model.getGravityVector()); %@@ requires theta=0
      gdz = g * s.Grad(Gt.cells.z);
   end

% --------------------------------------------------------------------%
   % function g = getGravityVector(model)
   %     % Get the gravity vector used to instantiate the model
   %     g = model.gravity;
   % end
end
end   
