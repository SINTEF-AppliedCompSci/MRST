function [Sw, krw, kro] = kr_correction(Sw, krw, kro)
%
% DESCRIPTION: checks the relative permeability table before before 
%              simulation - for things like sorting and if it goes from
%              zero to one
%
% SYNOPSIS:
%   [Sw, krw, kro] = kr_correction(Sw, krw, kro)
%
% PARAMETERS:
%   Sw: array water saturation
%   krw: array water relative permeability
%   kro: array oil relative permeability
%
% RETURNS:
%   Sw: array corrected water saturation
%   krw: array corrected water relative permeability
%   kro: array corrected oil relative permeability
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

% adding tale ends
if not(Sw(1) == 0)
    Sw = sort(Sw); krw = sort(krw);
    kro = sort(kro,'descend');
    Sw = [0;Sw]; krw = [0;krw];
    kro = [max(kro);kro];
end
if not(Sw(end) == 1)
    Sw = sort(Sw); krw = sort(krw);
    kro = sort(kro,'descend');
    Sw = [Sw;1]; krw = [krw;max(krw)];
    kro = [kro;0];
end

% check for invalid entries
krw(or(Sw>1,Sw<0)) = [];
kro(or(Sw>1,Sw<0)) = [];
Sw(or(Sw>1,Sw<0)) = [];
