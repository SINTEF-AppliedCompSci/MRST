function bG = boCO2(T_ref, rhoG, varargin)
% CO2 formation volume factor function
%
% SYNOPSIS:
%   function bG = boCO2(T_ref, rhoG, varargin)
%
% DESCRIPTION:
%
%   Produces a function representing the CO2 formation volume factor, which
%   expresses the density of CO2 as a fraction of some reference (surface)
%   density.  The returned function depends on pressure only (reference
%   temperature to be used is provided in the function call.
%
% PARAMETERS:
%   T_ref    - Reference temperature to be used (since the returned function
%              should only depend on pressure).  Since the reference
%              temperature may vary across the grid (due to a spatially
%              dependent temperature field, 'T_ref' can be given as a vector
%              of values.
%
%   rhoG     - Reference density
%
%   varargin - Optional parameters can be specified as key/value pairs on the
%              form ('key'/value ...).  These include:
%
%              * rho_datafile - name of datafile containing the sampled table to
%                               be used for computing CO2 densities as functions
%                               of pressure and temperature.
%
%             * sharp_phase_boundary - If 'true', will use one-sided evaluation
%                                      of derivatives near the liquid-vapor
%                                      boundary, in order to avoid smearing or
%                                      the derivatives across this
%                                      discontinuity.  Useful for plotting, but
%                                      in general not recommended for simulation
%                                      code using automatic differentiation,
%                                      since the discontinuous derivatives may
%                                      prevent the nonlinear solver from
%                                      converging.
%
% RETURNS:
%   bG - CO2 formation volume factor (a function of pressure)
%
% SEE ALSO:
%  `CO2props`, `SampledProp2D`

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

    opt.rho_datafile = 'CarbonDioxide_100000_400000000_278_524_800_800_D';
    opt.sharp_phase_boundary = true;
    opt = merge_options(opt, varargin{:});

    obj= CO2props('rhofile', opt.rho_datafile, ...
                  'sharp_phase_boundary', opt.sharp_phase_boundary);
    rhoCO2 =@(p) obj.rho(p, T_ref);
    bG =@(p, varargin) rhoCO2(p)/rhoG;
end
