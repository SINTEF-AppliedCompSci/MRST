function obj = CO2props(varargin)
% Generate a set of CO2 property functions based on sampled data.
%
% SYNOPSIS:
%   function obj = CO2props(varargin)
%
% DESCRIPTION:
%
%   Generates a structure containing a number of member functions representing
%   CO2 density, viscosity and enthalpy in terms of pressure and
%   temperature.  Derivative functions (and optionally, higher derivatives)
%   are also provided.  The functions are computed using sampled tables
%   (which can be changed), and support automatic differentiation as provided
%   within the MRST framework.
%
% PARAMETERS:
%   There are no required parameters.  A number of optional parameters can be
%   specified as key/value pairs on the form ('key'/value ...).  These
%   include:
%
%   - rhofile - location of sampled table with density values.  The default
%               table covers the pressure/temperature range [0.1, 400] MPa,
%               [278, 524] K. (400 x 400 samples).
%
%   - mufile  - location of sampled table with viscosity values.  The default
%               table covers the pressure/temperature range [0.1, 400] MPa,
%               [278, 524] K. (400 x 400 samples).
%
%   - hfile   - location of sampled table with enthalpy values.  The default
%               table covers the pressure/temperature range [0.1, 400] MPa,
%               [278, 524] K. (400 x 400 samples).
%
%   - const_derivatives - true or false.  If 'true', only first-order
%                         derivative functions will be included.
%
%   - assert  - if 'true', a property function will throw an error if user
%               tries to extrapolate outside the pressure/temperature range
%               covered by the corresponding sampled table.  If 'false',
%               either NaN values or extrapolated values will be returned in
%               this case, depending on whether the optional parameter
%               'nan_outside_range' is set to 'true' or 'false'.
%
%   - nan_outside_range - See documentation of 'assert' option above.
%
%   - sharp_phase_boundary - If 'true', will use one-sided evaluation of
%                            derivatives near the liquid-vapor boundary, in
%                            order to avoid smearing or the derivatives across
%                            this discontinuity.  Useful for plotting, but in
%                            general not recommended for simulation code
%                            using automatic differentiation, since the
%                            discontinuous derivatives may prevent the
%                            nonlinear solver from converging.
%
% RETURNS:
%   obj - Object containing a full set of CO2 property functions.
%
% SEE ALSO:
%   `SampledProp2D`

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

%% 
   opt.rhofile              = 'sampled_tables/CarbonDioxide_100000_400000000_278_524_800_800_D.mat';
   opt.mufile               = 'sampled_tables/CarbonDioxide_100000_400000000_278_524_800_800_V.mat';
   opt.hfile                = 'sampled_tables/CarbonDioxide_100000_400000000_278_524_800_800_H.mat';
   opt.const_derivatives    = false;
   opt.assert               = true;
   opt.nan_outside_range    = true;
   opt.sharp_phase_boundary = true;
   opt = merge_options(opt, varargin{:});

   boundary = [];
   if opt.sharp_phase_boundary
      [p_c, t_c]= CO2CriticalPoint();
      boundary = {[p_c, t_c], @CO2VaporPressure};
   end

   setupfun = @(name, file) ...
       SampledProp2D(name, file, ...
                     'phase_boundary'   , boundary             , ...
                     'nan_outside_range', opt.nan_outside_range, ...
                     'const_derivatives', opt.const_derivatives, ...
                     'assert_in_range'  , opt.assert);

   obj = setupfun('rho', opt.rhofile);
   obj = add_fields(obj, setupfun('mu' , opt.mufile));
   obj = add_fields(obj, setupfun('h'  , opt.hfile));


   %%%%
   obj.beta  = @(P, T)   obj.rhoDP(P, T) ./ obj.rho(P, T);  % Compressibility coef.
   obj.gamma = @(P, T) - obj.rhoDT(P, T) ./ obj.rho(P, T);  % Coef. of thermal expansion

  % compressibility factor cross derivative (equal for both Pbeta and Tbeta)
  obj.betaDPDT = @(P,T) (obj.rhoDPT(P, T) + obj.rhoDP(P, T) .* obj.rhoDT(P, T)) ./ obj.rho(P, T);

end

% ----------------------------------------------------------------------------

function obj = add_fields(obj, other)

   to_add = fieldnames(other);
   for i = 1:numel(to_add)
      obj.(to_add{i}) = other.(to_add{i});
   end
end
