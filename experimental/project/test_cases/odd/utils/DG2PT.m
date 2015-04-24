function [P, T] = DG2PT(depth, tgrad, t0, seafloor_depth, include_atm, rhoBrine)
%
% SYNOPSIS:
%   function [P, T] = DG2PT(depth, tgrad, t0, seafloor_depth, include_atm)
%
% PARAMETERS:
%   depth          - meters below sea level (not sea floor)
%   tgrad          - temperature gradients, to be given in degrees per km
%                    (typically within the range of 15 deg/km to 60 deg/km)
%   t0             - temperature at sea floor (not sea level) (277.15°K
%                    equals 4°C)
%   seafloor_depth - depth of sea floor
%   include_atm    - 'true' if atmospheric pressure should be added to pressure
%   rhoBrine       - density of brine
%
% RETURNS:
%   P - pressure in the reservoir at the given depth
%   T - temperatur in the reservoir at the given depth
%
% include_atm: set to 'true' if atmospheric pressure should be added to pressure 

% in general:
%     t0:      should be 4 degrees (centigrade) or 277.15 (Kelvin)
%     tgrad:   typical range is 15 deg/km to 50 deg/km

% calculating hydrostatic pressure assuming constant water density
P = depth * rhoBrine * 9.8; 

if include_atm
    P = P + atm;
end

% calculating temperature at the corresponding depth
T = t0 + ((depth - seafloor_depth)/1000) .* tgrad; 



