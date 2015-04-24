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
%               * props -  3-component logical vector (true/false) indicating
%                          which property functions to include, on the form:
%                          [include_density, include_visc., include_enthalpy]
%
%               * fixedT - If empty (default), the returned properties will be
%                          functions of pressure AND temperature.  If
%                          'fixedT' is a value (or vector of values), the
%                          returned properties will be functions of pressure
%                          only, with temperature consided constant and equal
%                          to the value(s) provided in 'fixedT'.
%
%                - 'assert_in_range':   If 'true', throw an error if user tried
%                                       to extrapolate outside valid range.
%                                       If 'false' (default), behavior in
%                                       this case will depend on the value of
%                                       'nan_outside_range'.
%
%                - 'nan_outside_range': If 'true' (default), return NaN
%                                       values outside valid range.
%                                       Otherwise, extrapolate as constant
%                                       function.
%
% RETURNS:
%   fluid - Fluid object endowed with the property functions (or a subset
%   thereof): rhoG, rhoW, muG, muW, hG, hW.
%
% SEE ALSO:
%  makeVEFluid, SampledProp2D, CO2props


   opt.pspan  = [0.1, 400] * mega * Pascal; % CO2 default pressure range
   opt.tspan  = [  4, 250] + 274;           % CO2 default temperature range
   opt.pnum   = 800; % number of pressure samples
   opt.tnum   = 800; % number of temperature samples
   opt.props  = [true false false]; % which props to include [rho, mu, h]
   opt.fixedT = [];
   opt.assert_range = false;
   opt.nan_outside_range = false;

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
      fluid.(['rho', shortname]) = load_property(opt, 'D', fluidname, opt.fixedT, ...
                                                 opt.assert_range, opt.nan_outside_range);
   end
   if opt.props(2)
      fluid.(['mu' , shortname]) = load_property(opt, 'V', fluidname, opt.fixedT, ...
                                                 opt.assert_range, opt.nan_outside_range);
   end
   if opt.props(3)
      fluid.(['h'  , shortname]) = load_property(opt, 'H', fluidname, opt.fixedT, ...
                                                 opt.assert_range, opt.nan_outside_range);

      if opt.props(1) % we have both enthalpy and density - we can also
                      % include internal energy
         fluid.(['u', shortname]) = ...
             @(P, T) fluid.(['h', shortname])(P, T) - P./fluid.(['rho',shortname])(P, T);
      end
   end
end

% ----------------------------------------------------------------------------

function pfun = load_property(opt, pname, fluidname, fixedT, assert_range, nan_outside)

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
                       'assert_in_range', assert_range, ...
                       'nan_outside_range', nan_outside);

   % We return the main evaluator function (which also works in an ADI-setting)
   pfun = obj.([pname]);%#ok

   if ~isempty(fixedT)
      % Temperature should be considered fixed -> property becomes function
      % of pressure only.
      pfun = @(p) pfun(p, fixedT);
   end

end
