function [Sw_pc, pc] = modified_skjaeveland_pc(swc, sor, cwi, coi, ri, ...
    bi, swd, sod)

% DESCRIPTION: create modified skjaeveland capillary pressure table
%
% SYNOPSIS:
%   [sw_pc, pc] = modified_skjaeveland_pc(swc, sor, cwi, coi, ri, bi,
%       swd, sod)
%
% PARAMETERS:
%   swc: connate water saturation
%   sor: residual oil saturation
%   cwi: water entry pressure (bar)
%   coi: oil entry pressure (bar)
%   ri: linear domain slope (bar/-)
%   bi: linear domain offset (bar)
%   swd: water domain saturation
%   sod: oil domain saturation
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

cwi = cwi * Convert('bar'); coi = coi * Convert('bar'); 
ri = ri * Convert('bar'); bi = bi * Convert('bar'); 

assert(swd > swc,'Swd must be higher than Swc!');
assert(sod < 1-sor, 'Sod must be lower than 1-Sor!');
Sw_pc = round(linspace(swc, 1 - sor, nPoints)',6);
pc = zeros(numel(Sw_pc),1);
Sw_pc(1) = Sw_pc(1) + eps; Sw_pc(end) = Sw_pc(end) - eps;
pc(Sw_pc <= swd) = cwi ./ ((Sw_pc(Sw_pc <= swd)-swc)./(swd-swc))...
            .^2 - cwi + swd * ri + bi;
pc(Sw_pc > swd & Sw_pc < sod)= Sw_pc(Sw_pc > swd & Sw_pc < sod) * ri + bi;  
pc(Sw_pc >= sod) = coi ./ ((1 - Sw_pc(Sw_pc >= sod) - sor) ./...
        (1 - sod - sor)) .^ 2 - coi + sod * ri + bi;
