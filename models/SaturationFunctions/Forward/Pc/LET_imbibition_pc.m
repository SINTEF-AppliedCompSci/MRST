function [sw_pc, pc] = LET_imbibition_pc(swc, sor, spontaneous_multiplier,...
    forced_multiplier, min_pc, max_pc, L_spont, E_spont, T_spont, L_forced, ...
    E_forced, T_forced, sw_pc0)
%
% DESCRIPTION: create LET imbibition capillary pressure table
%
% SYNOPSIS:
%   [Sw_pc, pc] = LET_imbibition_pc(swc, sor, spontaneous_multiplier,...
%     forced_multiplier, min_pc, max_pc, L_spont, E_spont, T_spont, L_forced, ...
%     E_forced, T_forced, sw_pc0)
%
% PARAMETERS:
%   swc: connate water saturation
%   sor: residual oil saturation
%   spontaneous_multiplier: activation for the spontaneous part (0 or 1)
%   forced_multiplier: activation for the forced part (0 or 1)
%   min_pc: minimum capillary pressure (bar)
%   max_pc: maximum capillary pressure (bar)
%   L_spont: spontaneous capillary pressur L exponent
%   E_spont: spontaneous capillary pressur E exponent
%   T_spont: spontaneous capillary pressur T exponent
%   L_forced: forced capillary pressur L exponent
%   E_forced: forced capillary pressur E exponent
%   T_forced: forced capillary pressur T exponent
%   sw_pc0: water saturation at capillary pressure of 0
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

min_pc = min_pc * Convert('bar');
max_pc = max_pc * Convert('bar');

Sw = linspace(0,1,nPoints)';
pc_spont = (spontaneous_multiplier .* max_pc .* 1e5 .* (1 - Sw) .^ L_spont) ./ ...
    ((1 - Sw) .^ L_spont + E_spont .* Sw .^ T_spont);
Se_spont = Sw * (sw_pc0 - swc) + swc;
pc_forced = (forced_multiplier .* min_pc .* 1e5 .* Sw .^ L_forced) ./ ...
    (Sw .^ L_forced + E_forced * (1 - Sw) .^ T_forced);
Se_forced = Sw * ((1 - sor) - sw_pc0) + sw_pc0;
if and(spontaneous_multiplier,forced_multiplier)
    sw_pc = [Se_spont(1:end-1); Se_forced];
    pc = [pc_spont(1:end-1); pc_forced]; 
elseif spontaneous_multiplier && not(forced_multiplier)
    sw_pc = Se_spont;
    pc = pc_spont; 
elseif forced_multiplier && not(spontaneous_multiplier)
    sw_pc = [Se_spont(1); Se_forced];
    pc = [pc_spont(1); pc_forced];         
end
