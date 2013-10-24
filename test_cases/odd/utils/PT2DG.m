function [D, G] = PT2DG(P, T, T0, seafloor_depth, include_atm)
%
% SYNOPSIS:
%   function [D, G] = PT2DG(P, T, T0, seafloor_depth, include_atm)
%
% PARAMETERS:
%   P              - Pressure value
%   T              - Temperature value
%   T0             - temperature at sea floor, in Kelvin (277.15°K equals
%                    4°C)
%   seafloor_depth - sea floor depth
%   include_atm    - does the value for pressure 'P' include the atmospheric
%                    pressure?
%
% RETURNS:
%   D - depth below sea level (not sea floor)
%   G - temperature gradient (degrees per kilometer)
%
% EXAMPLE:
%
% SEE ALSO:
%


% Convert from pressure/temperature to depth/temperature gradient

% depth: should be given in "meters below sea floor"
% tgrad: should be given in "degrees (of change) per 1000 m"

% in general:
%     t0:  should be 4 degrees (centigrade) or 277.15 (Kelvin)
%     tgrad: typical range is 15 deg/km to 50 deg/km

if include_atm
    P = P - atm;
end

% calculating depth
D = P / (9.8 * 1000);

% calculating temperature gradient
G = (T(:) - T0)./((D(:) - seafloor_depth)/1000);
