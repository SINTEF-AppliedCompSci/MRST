function [sw_kr, krw, kro] = burdine_kr(sw_pc, pc, krwSor, kroSwc)
% <keywords>
%
% Purpose : create a relative permeability table based on burdine integral
%
% Syntax :
%   [Sw_pc, krw, kro] = burdine_kr(Sw_pc, pc, krwSor, kroSwc)
%
% Input Parameters :
%   Sw_pc: water saturation array
%   pc: capillary pressure array corresponding to Sw_pc array (in pascal)
%   krwSor: end point water relative permeability
%   kroSwc: end point oil relative permeability
%
% Return Parameters :
%   Sw_kr: water saturation array for the relative permeability table
%           this is the same as the input water saturation array
%   krw: water relative permeability array
%   kro: oil relative permeability array
%
% Description :
%
% Author : 
%    Siroos Azizmohammadi
%    Omidreza Amrollahinasab
%
% History :
% \change{1.0}{23-Nov-2021}{Original}
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
end
