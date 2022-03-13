classdef CO2VEBlackOilTypeModelSens < ReservoirModel
    % Black-oil type model for vertically integrated gas/water flow
    % 
    % SYNOPSIS:
    %    model = CO2VEBlackOilTypeModel(Gt, rock2D, fluid, varargin)
    %
    % DESCRIPTION:
    %   Class representing a model with vertically-integrated two-phase flow (CO2
    %   and brine), based on the s-formulation where upscaled saturation is a
    %   primary variable), and with optional support for dissolution of gas into
    %   the water phase.
    %
    % REQUIRED PARAMETERS:
    %   Gt     - Top surface grid, generated from a regular 3D simulation grid
    %            using the 'topSurfaceGrid' function in MRST-co2lab
    %   
    %   rock2D - Vertically averaged rock structure, generated from regular 3D
    %            rock structure using the 'averageRock' function in MRST-co2lab.
    % 
    %   fluid  - Fluid object, representing the properties of the water and gas
    %            phases.  The fluid object can be constructed using the
    %            'makeVEFluid' function in MRST-co2lab.  This object also
    %            specifies whether or not gas dissolves into water, and if so,
    %            whether to model dissolution as an instantaneous or rate-driven
    %            process. 
    % RETURNS:
    %   Class instance.
    %
    % SEE ALSO:
    % `topSurfaceGrid`, `averageRock`, `makeVEFluid`, `ReservoirModel`

    % ============================= Class properties ==========================
properties
   
   % Equation is chosen based on whether fluid object includes dissolution
   % effects or not 
   equation
   
end
      
% ============================== Public methods ===========================
methods
   
   % Constructor
   function model = CO2VEBlackOilTypeModelSens(Gt, rock2D, fluid, varargin)
   
      opt = struct();
      [opt, unparsed] = merge_options(opt, varargin{:});%#ok
   
      model@ReservoirModel(Gt, varargin{:});
      model.rock    = rock2D;
      model.fluid   = fluid;
      model.water   = true;
      model.gas     = true;
      model.oil     = false;
      model.gravity = gravity;
      
      if isfield(fluid, 'dis_rate')
         % use model equations with dissolution
         model.equation = @equationsWGVEdisgas;
         error('sensitivites??')
      else
         % use basic model equations (no dissolution)
         model.equation = @equationsWGVEbasicSens;
         %model.adjointType=opt.adjointType;
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

      % Compute vertially-integrated transmissibilities if not provided
      rock_tmp      = rock; 
      rock_tmp.perm = rock.perm .* Gt.cells.H; 
      T             = computeTrans(Gt, rock_tmp); 
      cf            = Gt.cells.faces(:, 1); 
      nf            = Gt.faces.num; 
      T             = 1 ./ accumarray(cf, 1 ./ T, [nf, 1]);       

      % Computing vertically-integrated pore - volume
      pv = poreVolume(Gt, rock); 
      pv = pv .* Gt.cells.H; 

      model.operators = setupOperatorsTPFA(Gt, rock, 'porv', pv, 'trans', T);
      model.operators.T_all = T;
      
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
function forces = getValidDrivingForces(model) 
        forces = struct('W', [], 'bc', [], 'src', [],'dz',[],'rhofac',1,'permfac',1,'porofac',1);
end

% ------------------------------------------------------------------------
   function [fn, index] = getVariableField(model, name, varargin)
      
      switch(lower(name))
        case {'sgmax'}
          index = 1;
          fn = 'sGmax';
        case {'rs'}
          index = 1;
          fn = 'rs';
        case {'dz'}
            index = 1;
            fn = 'dz';
        case {'rhofac'}
            index = 1;
            fn = 'rhofac';
        case {'permfac'}
            index = 1;
            fn = 'permfac';
        case {'porofac'}
            index = 1;
            fn = 'porofac';
        otherwise
          [fn, index] = getVariableField@ReservoirModel(model, name, varargin{:});
      end
   end

% ----------------------------------------------------------------------------
function [state, report] = updateState(model, state, problem, dx, drivingForces)

   [state, report] = updateState@ReservoirModel(model, state, problem, dx, ...
                                                drivingForces);
   %probably should limit update of dz
                                            
   if isfield(model.fluid, 'dis_rate')
       error('')
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
      g  = model.gravity(3); % @@ requires theta=0
      gdz = g * s.Grad(Gt.cells.z);
   end

end
end   
%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

