function [sw_kr, krw, kro] = van_genuchten_burdine_kr(swc, sor, krwSor, ...
    kroSwc, n_kr)
%
% DESCRIPTION: create a relative permeability table based on van
%              genuchten burdine integral
%
% SYNOPSIS:
%   van_genuchten_burdine_kr(swc, sor, krwSor, ...
%    kroSwc, n_kr)
%
% PARAMETERS:
%   swc: connate water saturation
%   sor: residual oil saturation
%   krwSor: end point water relative permeability
%   kroSwc: end point oil relative permeability
%   n_kr: van genuchten exponent
%
% RETURNS:
%   Sw_kr: water saturation array for the relative permeability table
%           this is the same as the input water saturation array
%   krw: water relative permeability array
%   kro: oil relative permeability array
%
% ----------------------------------
% (c) 2020-2022
% Siroos Azizmohammadi
% Omidreza Amrollahinasab
% Montanuniversität Leoben, Austria
% Chair of Reservoir Engineering
% https://dpe.ac.at/
% ----------------------------------
%
%%
assert(swc >= 0 & swc < 1,"Wrong input for Swc")
assert(sor >= 0 & sor < 1-swc,"Wrong input for Sor")
nPoints = 50;
sw  = round(linspace(swc, 1 - sor, nPoints)',6);
sw_star  = round((sw - swc) ./ (1 - swc - sor),6);
sw_star(1) = 0; sw_star(end) = 1;

m = (n_kr - 1) / n_kr;
krw = krwSor * (sw_star .^ 2) .* (1 - (1 - sw_star .^ (1 / m)) .^ m);
kro = kroSwc * ((1 - sw_star) .^ 2) .* (1 - sw_star .^ (1 / m)) .^ m;
[sw_kr, krw, kro] = kr_correction(sw, krw, kro);
