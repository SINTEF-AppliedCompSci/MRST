function [sw_pc, pc] = modified_skjaeveland_masalmeh_pc(swc, sor, cwi, coi, awi, aoi, sw_cutoff, b)
% <keywords>
%
% Purpose : create modified skjaeveland from masalmeh capillary pressure table
%
% Syntax :
%   [Sw_pc, pc] = modified_skjaeveland_masalmeh_pc(Swc, Sor, cwi, coi, awi, aoi, sw_cutoff, b)
%
% Input Parameters :
%   Swc: connate water saturation
%   Sor: residual oil saturation
%   cwi: water entry pressure (bar)
%   coi: oil entry pressure (bar)
%   awi: water exponent
%   aoi: oil exponent
%   sw_cutoff: sw cutoff
%   b: cutoff multiplier (bar)
%
% Return Parameters :
%   Sw_pc: array water saturation
%   pc: array capillary pressure (pascal)
%
% Description :
%
% Author : 
%    Siroos Azizmohammadi
%    Omidreza Amrollahinasab
%
% History :
% \change{1.0}{09-Nov-2021}{Original}
%
% --------------------------------------------------
% (c) 2021, Siroos Azizmohammadi,
% Omidreza Amrollahinasab
% Chair of Reservoir Engineering, University of Leoben, Austria
% email: info@dpe.ac.at
% url: dpe.ac.at
% --------------------------------------------------
%
%%
    assert(swc >= 0 & swc < 1,"Wrong input for Swc")
    assert(sor >= 0 & sor < 1-swc,"Wrong input for Sor")
    nPoints = 50;
    eps = 1e-2; % To avoind infinite values
    
    cwi = cwi * Convert('bar');
    coi = coi * Convert('bar');
    b = b * Convert('bar');
    
    sw_pc = round(linspace(swc, 1 - sor, nPoints)',6);
    Sew  = round((sw_pc - swc) ./ (1 - swc),6);
    Sew(1) = eps; Sew(end) = 1; So_pc = 1 - sw_pc;
    Seo  = round((So_pc - sor) ./ (1 - sor),6);
    Seo(1) = 1; Seo(end) = eps;
    pc = cwi ./ (Sew .^ awi) + coi ./ (Seo .^ aoi) + b * (sw_cutoff - sw_pc);
end
