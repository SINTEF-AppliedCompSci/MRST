function [theta, krw, C_theta] = vGM_theta(phys)
% van Genuchten-Mualem retention curves for psi-theta model.
%
% SYNOPSIS:
%   [theta, krw, C_theta] =  vGM_theta(phys)
%
% PARAMETERS:
%   phys     - Structure, containing physical parameters
%
%  RETURNS:
%   theta    - Function, Water content, theta = theta(psi)
%   krw      - Function, Relative permeability, krw = krw(psi)
%   C_theta  - Function, Specific moisture capacity, C_theta = C_theta(psi)
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
alpha = phys.flow.alpha;
theta_r = phys.flow.theta_r;
theta_s = phys.flow.theta_s;
n = phys.flow.n;
m = phys.flow.m;

% Boolean inline function: true for unsat zone, false for saturated zone
isUnsat = @(psi) psi < 0;

% Water content
theta= @(psi)   isUnsat(psi) .* ((theta_s - theta_r) ...
                ./ (1 + (alpha .* abs(psi)) .^ n) .^ m + theta_r ) ...
                + ~isUnsat(psi) .* theta_s;
    
% Specific Moisture Capacity
C_theta= @(psi)   isUnsat(psi) .* ((m .* n .* psi .* (theta_r - theta_s)  ...
                .* alpha .^ n .* abs(psi) .^ (n - 2)) ...
                ./ (alpha ^ n .* abs(psi) .^ n + 1) .^ (m+1)) ...
                + ~isUnsat(psi) .* 0;

% Relative permeability
krw= @(psi)   isUnsat(psi) .* ((1 - (alpha .* abs(psi)) .^ (n - 1) ...
              .* (1 + (alpha .* abs(psi)) .^ n) .^ (-m)) .^ 2 ...
              ./ (1 + (alpha .* abs(psi)) .^ n) .^ (m ./ 2)) ...
              + ~isUnsat(psi) .* 1;