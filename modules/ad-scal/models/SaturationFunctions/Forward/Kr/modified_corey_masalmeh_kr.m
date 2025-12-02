function [sw, krw, kro] = modified_corey_masalmeh_kr(swc, sor, krwSor, ...
    kroSwc, nW, nNW, cW, cNW)
%
% DESCRIPTION: create modified corey masalmeh relative permeability table
%
% SYNOPSIS:
%   [sw, krw, kro] = modified_corey_masalmeh_kr(swc, sor, krwSor, ...
%   kroSwc, nW, nNW, cW, cNW)
%
% PARAMETERS:
%   Swc: connate water saturation
%   Sor: residual oil saturation
%   krwSor: end point water relative permeability
%   kroSwc: end point oil relative permeability
%   nW: wetting phase corey exponent
%   nNW: non-wetting phase corey exponent
%   cW: water tail exponent
%   cNW: non-wetting tail exponent
%
% RETURNS:
%   Sw: array water saturation
%   krw: array water relative permeability
%   kro: array oil relative permeability
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
Se  = round((sw - swc) ./ (1 - swc - sor),6);
Se(1) = 0; Se(end) = 1;
krw = (krwSor .* (Se) .^ nW) + (cW / (1+cW) .* Se);
kro = (kroSwc .* (1 - Se) .^ nNW) + (cNW /...
    (1+cNW) .* (1 - Se));
krw(1) = 0; kro(end) = 0;
[sw, krw, kro] = kr_correction(sw, krw, kro);
