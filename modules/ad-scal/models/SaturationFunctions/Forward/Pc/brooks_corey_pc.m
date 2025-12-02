function [sw_pc, pc] = brooks_corey_pc(swc_pc, sor_pc, pd, lambda)
%
% DESCRIPTION: create a capillary pressure table using brooks corey model
%
% SYNOPSIS:
%   [sw_pc, pc] = brooks_corey_pc(swc_pc, sor_pc, pd, lambda)
%
% PARAMETERS:
%   swc: connate water saturation
%   sor: residual oil saturation
%   pd: entry pressure (bar)
%   lambda: exponent
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
assert(swc_pc >= 0 & swc_pc < 1,"Wrong input for Swc")
assert(sor_pc >= 0 & sor_pc < 1-swc_pc,"Wrong input for Sor")
nPoints = 50;
eps = 1e-2; % To avoind infinite values


pd = pd * Convert('bar');
sw_pc  = round(linspace(swc_pc, 1 - sor_pc, nPoints)',6);
Se  = round((sw_pc - swc_pc) ./ (1 - swc_pc - sor_pc),6);
Se(1) = eps; Se(end) = 1;
pc = pd * (Se) .^ (-1 / lambda);
