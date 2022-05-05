function [sw_kr, krw, kro] = burdine_kr(sw_pc, pc, krwSor, kroSwc)
%
% DESCRIPTION: create a relative permeability table based on burdine integral
%
% SYNOPSIS:
%   [Sw_pc, krw, kro] = burdine_kr(Sw_pc, pc, krwSor, kroSwc)
%
% PARAMETERS:
%   Sw_pc: water saturation array
%   pc: capillary pressure array corresponding to Sw_pc array (in pascal)
%   krwSor: end point water relative permeability
%   kroSwc: end point oil relative permeability
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
%avoiding problems with the integral 
pc(pc==0) = min(pc(pc>0));

Swc = min(sw_pc); Sor = 1 - max(sw_pc);
sw_star = round((sw_pc - Swc)./(1- Swc - Sor),6);
burdine_pc = 1./pc.^2;
F = griddedInterpolant(sw_star,burdine_pc);
fun = @(t) F(t);
q = integral(fun, sw_star(1), sw_star(end));
krw = []; kro = [];
for i = 1: length(sw_star)
    krw = [krw; krwSor .* (sw_star(i)).^2 .* ...
        integral(fun, sw_star(1), sw_star(i)) ./ q];
    kro = [kro; kroSwc .* (1 - sw_star(i)).^2 .* ...
        integral(fun, sw_star(i), sw_star(end)) ./ q];
end
[sw_kr, krw, kro] = kr_correction(sw_pc, krw, kro);

