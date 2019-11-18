function [S_w, krw, C_S] = vanGenuchtenMualemSw(a, S_r, n, m)
% van Genuchten-Mualem retention curves for pw-Sw model.
%
% SYNOPSIS:
%   [S_w, krw, C_S] = vanGenuchtenMualemSw(a, S_r, n, m)
%
% PARAMETERS:
%   a      - Scalar, a = alphaVan/gamma (van Genuchten equation parameter)
%   S_r    - Scalar, residual saturation
%   n      - Scalar, van Genuchten equation parameter
%   m      - Scalar, van Genuchten equation parameter
%
%  RETURNS:
%   S_w    - Function, Saturation, i.e., S_w = S_w(p_w)
%   krw    - Function, Relative permeability, i.e., krw = krw(p_w)
%   C_S    - Function, Specific saturation capacity, i.e., C_S = C_S(p_w)
%

%{
Copyright 2018-2019, University of Bergen.

This file is part of the fv-unsat module.

fv-unsat is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

fv-unsat is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%} 

% Boolean function that determines if we are in the 
% unsat (true) or sat (false) zone
isUnsat = @(p) p < 0;

% Saturation
S_w = @(p)  isUnsat(p) .* (((1-S_r)./(1+(a.*abs(p)).^n).^m)+S_r) ...
            + ~isUnsat(p) .* 1;

% Specific saturation capacity
C_S = @(p)  isUnsat(p) .* ((-m.*n.*(1-S_r).*p.*(a.*abs(p)).^n) ...
            ./(abs(p).^2.*(1+(a.*abs(p)).^n).^(m+1))) ... 
            + ~isUnsat(p) .* 0;

% Relative permeability
krw = @(p)  isUnsat(p) .* ((1-(a.*abs(p)).^(n-1) ...
            .*(1+(a.*abs(p)).^n).^(-m)).^2 ...
            ./(1+(a.*abs(p)).^n).^(m./2)) ...
            + ~isUnsat(p) .* 1;