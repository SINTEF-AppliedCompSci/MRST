function [theta, krw, C_theta] = vanGenuchtenMualemTheta(alpha, theta_s, theta_r, n, m)
% van Genuchten-Mualem retention curves for psi-theta model.
%
% SYNOPSIS:
%   [theta, krw, C_theta] = vanGenuchtenMualemTheta(alpha, theta_s, theta_r, n, m)
%
% PARAMETERS:
%   alpha    - Scalar, van Genuchten equation parameter
%   theta_s  - Scalar, water content at saturated conditions
%   theta_r  - Scalar, residual water content
%   n        - Scalar, van Genuchten equation parameter
%   m        - Scalar, van Genuchten equation parameter
%
%  RETURNS:
%   theta    - Function, Water content, i.e., theta = theta(psi)
%   krw      - Function, Relative permeability, i.e., krw = krw(psi)
%   C_theta  - Function, Specific saturation capacity, i.e., C_theta = C_theta(psi)
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