function [sw, krw, kro] = LET_kr(swc, sor, krwSor, kroSwc, Lw, Ew, ...
    Tw, Lnw, Enw, Tnw)
%
% DESCRIPTION: create LET relative permeability table
%
% SYNOPSIS:
%   [sw, krw, kro] = LET_kr(swc, sor, krwSor, kroSwc, Lw, Ew, ...
%    Tw, Lnw, Enw, Tnw)
%
% PARAMETERS:
%   Swc: connate water saturation
%   Sor: residual oil saturation
%   krwSor: end point water relative permeability
%   kroSwc: end point oil relative permeability
%   Lw: wetting phase L exponent
%   Ew: wetting phase E exponent
%   Tw: wetting phase T exponent
%   Lnw: non-wetting phase L exponent
%   Enw: non-wetting phase E exponent
%   Tnw: non-wetting phase T exponent
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
krw = krwSor * ( (Se .^ Lw) ./ (Se .^ Lw + Ew .*...
    (1 - Se) .^ Tw));
kro = kroSwc * ( ( (1 - Se) .^ Lnw ) ./ ( (1 - Se) .^ Lnw...
    + Enw .* (Se .^ Tnw) ) );
krw(1) = 0; kro(end) = 0;
[sw, krw, kro] = kr_correction(sw, krw, kro);

