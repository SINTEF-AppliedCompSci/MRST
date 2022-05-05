function [sw_pc, pc] = skjaeveland_pc(swc, sor, cwi, coi, awi, aoi)
%
% DESCRIPTION: create skjaeveland capillary pressure table
%
% SYNOPSIS:
%   [sw_pc, pc] = skjaeveland_pc(swc, sor, cwi, coi, awi, aoi)
%
% PARAMETERS:
%   swc: connate water saturation
%   sor: residual oil saturation
%   cwi: water entry pressure (bar)
%   coi: oil entry pressure (bar)
%   awi: water exponent
%   aoi: oil exponent
%
% RETURNS:
%   sw_pc: array water saturation
%   pc: array capillary pressure (pascal)
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
eps = 1e-2; % To avoind infinite values

cwi = cwi * Convert('bar');
coi = coi * Convert('bar');
sw_pc = round(linspace(swc, 1 - sor, nPoints)',6);
Sew  = round((sw_pc - swc) ./ (1 - swc),6);
Sew(1) = eps; Sew(end) = 1; So_pc = 1 - sw_pc;
Seo  = round((So_pc - sor) ./ (1 - sor),6);
Seo(1) = 1; Seo(end) = eps;
pc = cwi ./ (Sew .^ awi) + coi ./ (Seo .^ aoi);
