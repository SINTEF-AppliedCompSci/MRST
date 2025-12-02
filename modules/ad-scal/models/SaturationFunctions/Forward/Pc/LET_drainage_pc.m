function [sw_pc, pc] = LET_drainage_pc(swc, entry_multiplier,...
    forced_multiplier, entry_pc, max_pc, L_entry, E_entry, T_entry, L_forced, ...
    E_forced, T_forced)
%
% DESCRIPTION: create LET drainage capillary pressure table
%
% SYNOPSIS:
%   [Sw_pc, pc] = LET_drainage_pc(swc, entry_multiplier,...
%     forced_multiplier, entry_pc, max_pc, L_entry, E_entry, T_entry, L_forced, ...
%     E_forced, T_forced)
%
% PARAMETERS:
%   swc: connate water saturation
%   sor: residual oil saturation
%   entry_multiplier: activation for the entry part (0 or 1)
%   forced_multiplier: activation for the forced part (0 or 1)
%   entry_pc: entry capillary pressure (bar)
%   max_pc: maximum capillary pressure (bar)
%   L_entry: entry capillary pressur L exponent
%   E_entry: entry capillary pressur E exponent
%   T_entry: entry capillary pressur T exponent
%   L_forced: forced capillary pressur L exponent
%   E_forced: forced capillary pressur E exponent
%   T_forced: forced capillary pressur T exponent
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
nPoints = 50;

Sw = linspace(0,1,nPoints)';
max_pc = max_pc * Convert('bar');
entry_pc = entry_pc * Convert('bar');
Se_entry = swc + Sw .* (1 - swc);
pc_entry = (entry_pc .* Sw .^ L_entry) ./ ...
    (Sw .^ L_entry + E_entry .* (1 - Sw) .^ T_entry);
pc_forced = forced_multiplier .* (max_pc - entry_pc) .* ...
    (((1- Sw) .^ T_forced) ./ ((1 - Sw) .^ T_forced + E_forced .* Sw .^ T_forced)) .^ (1 / L_forced) ...
    - entry_multiplier .* pc_entry + entry_pc;
sw_pc = Se_entry; pc = pc_forced;

