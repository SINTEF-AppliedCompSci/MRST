function fluid = addSampledFluidProperties(fluid, shortname, varargin)
% Add density, viscosity or enthalpy properties to a fluid object.  The
% properties may be those of CO2 or of water.
%
% SYNOPSIS:
%   function fluid = addSampledFluidProperties(fluid, shortname, varargin)
%
% DESCRIPTION:
%   This function endows a preexisting fluid object with property functions
%   for density, viscosity and/or enthalpy.  These functions will depend on
%   pressure and temperature, or on pressure only (if requested), and are
%   computed based on sampled tables.
%
% PARAMETERS:
%   fluid     - Preexisting fluid object to receive the requested property
%               functions.
%
%   shortname - Either 'G' (for gas, i.e. CO2 property functions) or 'W', for
%               water property functions.
%
%   varargin  - Optional arguments supplied as 'key'/value pairs ('pn'/nv).
%               These include:
%
%               * pspan -  2-component vector [pmin pmax] to specify the upper
%                          and lower pressure range (Pa) of the sampled table.
%                          If no such table preexists, it will be computed using
%                          the CoolProps software package (if this is not
%                          available, an error will be thrown).
%
%               * tspan -  2-component vector [tmin tmax] to specify the upper
%                          and lower temperature range (K) of the sampled
%                          table.  If no such table preexists, it will be
%                          computed using the CoolProps software package (if
%                          this is not available, an error will be thrown).
%
%               * pnum  -  number of (equidistant) samples along the pressure
%                          dimension in the sampled table.
%
%               * tnum  -  number of (equidistant) samples along the
%                          temperature dimension in the sampled table.
%
%               * props - 4-component logical vector (true/false) indicating
%                          which property functions to include, on the form:
%                          [density, viscosity, enthalpy, conductivity]
%
%               * fixedT - If empty (default), the returned properties will be
%                          functions of pressure AND temperature.  If
%                          'fixedT' is a value (or vector of values), the
%                          returned properties will be functions of pressure
%                          only, with temperature consided constant and equal
%                          to the value(s) provided in 'fixedT'.
%
%               * 'assert_in_range':   If 'true', throw an error if user tried
%                                      to extrapolate outside valid range.
%                                      If 'false' (default), behavior in
%                                      this case will depend on the value of
%                                      'nan_outside_range'.
%
%               * 'nan_outside_range': If 'true' (default), return NaN
%                                      values outside valid range.
%                                      Otherwise, extrapolate as constant
%                                      function.
%           
%               * 'include_derivatives': include explicit functions for
%                                        partial derivatives in pressure and 
%                                        temperature (default: false)
%
% RETURNS:
%   fluid - Fluid object endowed with the property functions (or a subset
%   thereof): rhoG, rhoW, muG, muW, hG, hW.
%
% SEE ALSO:
%  `makeVEFluid`, `SampledProp2D`, `CO2props`

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
   opt.pspan  = [0.1, 400] * mega * Pascal; % CO2 default pressure range
   opt.tspan  = [  4, 250] + 274;           % CO2 default temperature range
   opt.pnum   = 800; % number of pressure samples
   opt.tnum   = 800; % number of temperature samples
   opt.props  = [true false false false]; % which props to include [rho, mu, h, lambda]
   opt.fixedT = [];
   opt.assert_in_range = false;
   opt.nan_outside_range = false;
   opt.include_derivatives = false;

   opt = merge_options(opt, varargin{:});

   % Determine full fluid name
   switch shortname
     case 'G'
       fluidname = 'CarbonDioxide';
     case 'W'
       fluidname = 'Water';
     otherwise
      error('Unsupported fluid name');
   end

   % Add density, viscosity and enthalpy properties
   if opt.props(1)
      [f, fdp, fdt] = load_property(opt, 'D', fluidname, opt.fixedT, ...
                                    opt.assert_in_range, opt.nan_outside_range);
      fname = ['rho', shortname];
      
      fluid.(fname) = f;
      if opt.include_derivatives
         [fluid.([fname, '_dp']), fluid.([fname, '_dt'])] = deal(fdp, fdt);
      end
   end
   if opt.props(2)
      [f, fdp, fdt] = load_property(opt, 'V', fluidname, opt.fixedT, ...
                                    opt.assert_in_range, opt.nan_outside_range); 
      
      fname = ['mu', shortname];
      
      fluid.(fname) = f;
      if opt.include_derivatives
         [fluid.([fname, '_dp']), fluid.([fname, '_dt'])] = deal(fdp, fdt);
      end
                                                 
   end
   if opt.props(3)
      [f, fdp, fdt] = load_property(opt, 'H', fluidname, opt.fixedT, ...
                                    opt.assert_in_range, opt.nan_outside_range); 
      
      fname = ['h', shortname];
      
      fluid.(fname) = f;
      if opt.include_derivatives
         [fluid.([fname, '_dp']), fluid.([fname, '_dt'])] = deal(fdp, fdt);
      end

      if opt.props(1) % we have both enthalpy and density - we can also
                      % include internal energy
         fluid.(['u', shortname]) = ...
             @(P, T) fluid.(['h', shortname])(P, T) - P./fluid.(['rho',shortname])(P, T);
      end
   end
   if (numel(opt.props) > 3 && opt.props(4))
      [f, fdp, fdt] = load_property(opt, 'L', fluidname, opt.fixedT, ...
                                    opt.assert_in_range, opt.nan_outside_range); 
      fname = ['lambda', shortname];
      
      fluid.(fname) = f;
      if opt.include_derivatives
         [fluid.([fname, '_dp']), fluid.([fname, '_dt'])] = deal(fdp, fdt);
      end
   end
end

% ----------------------------------------------------------------------------

function [fun, fun_dp, fun_dt] = load_property(opt, pname, fluidname, fixedT, assert_in_range, nan_outside)

   tabledir = [fileparts(mfilename('fullpath')) '/sampled_tables/'];
   fname = [tabledir, propFilename(opt.pspan, opt.tspan, opt.pnum, opt.tnum, fluidname, pname)];

   if (exist(fname) ~= 2)%#ok
      % data table not yet generated.  The following command will generate
      % and save them. 
      try
         generatePropsTable(tabledir, fluidname, pname, opt.pspan, opt.tspan, opt.pnum, ...
                            opt.tnum);
      catch ME
         if strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
            error(['Failed to generate property tables, as CoolProps could not ' ...
                   'be found.  Make sure you have CoolProps installed with ' ...
                   'the Matlab wrapper, and that the directory of the ' ...
                   'corresponding ''PropsSI'' is in your Matlab search ' ...
                   'path.']);
         else
            error(['Generating new property table failed due to an unknown ' ...
                   'error.']);
         end
      end
   end

   % We here load the generated table, and construct an object with evaluator
   % functions
   obj = SampledProp2D(pname, fname, ...
                       'assert_in_range', assert_in_range, ...
                       'nan_outside_range', nan_outside);

   % We return the main evaluator function (which also works in an ADI-setting)
   fun = obj.([pname]);%#ok
   fun_dp = obj.([pname, 'DP']);
   fun_dt = obj.([pname, 'DT']);

   if ~isempty(fixedT)
      % Temperature should be considered fixed -> property becomes function
      % of pressure only.
      fun = @(p) fun(p, fixedT);
      fun_dp = @(p) fun_dp(p, fixedT);
      fun_dt = []; 
   end

end
