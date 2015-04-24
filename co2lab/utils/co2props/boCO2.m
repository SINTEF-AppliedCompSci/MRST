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
%  CO2props, SampledProp2D
   
    opt.rho_datafile = 'rho_big_trunc';
    opt.sharp_phase_boundary = true;
    opt = merge_options(opt, varargin{:});
    
    obj= CO2props('rhofile', opt.rho_datafile, ...
                  'sharp_phase_boundary', opt.sharp_phase_boundary);
    rhoCO2 =@(p) obj.rho(p, T_ref);
    bG =@(p) rhoCO2(p)/rhoG;
end