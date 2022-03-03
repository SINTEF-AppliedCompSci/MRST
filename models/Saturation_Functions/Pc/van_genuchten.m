function [sw_pc, pc_array] = van_genuchten(swc, sor, alpha_pc, n_pc)
% <keywords>
%
% Purpose : create van genuchten capillary pressure table
%
% Syntax :
%   [Sw_pc, pc_array] = van_genuchten(Swc, Sor, alpha, n)
%
% Input Parameters :
%   swc: connate water saturation
%   sor: residual oil saturation
%   alpha_pc: entry pressure (1/bar)
%   n_pc: van genuchten exponent
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
   
    alpha_pc = alpha_pc / Convert('bar');
    sw_pc  = round(linspace(swc, 1 - sor, nPoints)',6);
    Se  = round((sw_pc - swc) ./ (1 - swc - sor),6);
    Se(1) = eps; Se(end) = 1;
    m = (n_pc - 1) / n_pc;
    pc_array = 1 / alpha_pc * ( ( (Se) .^ (-1 / m) - 1 ) .^ (1 / n_pc) );
end