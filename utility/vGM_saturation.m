function [S_w, krw, C_S] = vGM_saturation(phys)
% van Genuchten-Mualem retention curves for pw-Sw model.
%
% SYNOPSIS:
%   [S_w, krw, C_S] = vGM_saturation(phys)
%
% PARAMETERS:
%   phys   - Structure, containing physical parameters
%
%  RETURNS:
%   S_w    - Function, Saturation, S_w = S_w(p_w)
%   krw    - Function, Relative permeability, krw = krw(p_w)
%   C_S    - Function, Specific saturation capacity, C_S = C_S(p_w)
%

%{
Copyright 2018-2020, University of Bergen.

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

% Retrieving physical properties
a = phys.flow.alpha / (phys.flow.rho * phys.flow.g); 
S_r = phys.flow.theta_r / phys.flow.poro;
n = phys.flow.n;
m = phys.flow.m;


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