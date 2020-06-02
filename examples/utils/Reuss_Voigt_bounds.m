function [x_r, x_v] = Reuss_Voigt_bounds(x_m, x_f, vol_m, vol_f)
%
% SYNOPSIS:
%   function [x_r, x_v] = Reuss_Voigt_bounds(x_m, x_f, vol_m, vol_f)
%
% DESCRIPTION:
%   Function to calculate the upscaled properties of a two phase material
%   using the Reuss and Voigt (lower and upper resp.) of the given material
%   properties x_i.
%
% PARAMETERS:
%   x_m     - mechanical property of matrix continuum
%   x_f     - mechanical property of fracture continuum
%   vol_m   - matrix volume fraction
%   vol_f   - fracture volume fraction
%
% RETURNS:
%   x_r    - Reuss bound of mechanical property x
%   x_v    - Reuss bound of mechanical property x
%
% SEE ALSO: HS_bound
%
x_r = 1/(vol_m/x_m + vol_f/x_f);
x_v = x_m*vol_m + x_f*vol_f;

end
