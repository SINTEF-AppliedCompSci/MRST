function [Sw_pc, pc] = brooks_corey_pc(Swc_pc, Sor_pc, pd, lambda)
% <keywords>
%
% Purpose : create brooks corey capillary pressure table
%
% Syntax :
%   [Sw_pc, pc] = brooks_corey_pc(Swc, Sor, pd, lambda)
%
% Input Parameters :
%   Swc: connate water saturation
%   Sor: residual oil saturation
%   pd: entry pressure (bar)
%   lambda: exponent
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
    assert(Swc_pc >= 0 & Swc_pc < 1,"Wrong input for Swc")
    assert(Sor_pc >= 0 & Sor_pc < 1-Swc_pc,"Wrong input for Sor")
    nPoints = 50;
    eps = 1e-2; % To avoind infinite values
    
    
    pd = pd * Convert('bar');
    Sw_pc  = round(linspace(Swc_pc, 1 - Sor_pc, nPoints)',6);
    Se  = round((Sw_pc - Swc_pc) ./ (1 - Swc_pc - Sor_pc),6);
    Se(1) = eps; Se(end) = 1;
    pc = pd * (Se) .^ (-1 / lambda);
end
