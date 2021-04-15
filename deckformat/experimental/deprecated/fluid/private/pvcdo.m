function [B, dB, R, dR, mu, P] = pvcdo(tab, p)
%Evaluate dead oil properties for oil with constant compressibility.
%
% SYNOPSIS:
%   [B, dB, R, dR, mu, P] = pvcdo(tab, p)
%
% PARAMETERS:
%
% RETURNS:
%   B, dB - Volume formation factor (B) and derivative of B with respect to
%           pressure (dB).
%
%   R, dR - Miscibility ratio (R) and derivative of R with respect to
%           pressure (dR).  This is a dead oil table so R == dR == 0.
%
%   mu    - Oil viscosity.
%
%   P     - Predicate for identifying which cells are in a saturated state.
%
% SEE ALSO:
%   `pvtg`, `pvto`.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}


   p0  = tab(1);   % Reference pressure.
   B0  = tab(2);   % Oil formation factor at reference pressure.
   C   = tab(3);   % Oil compressibility (constant).
   mu0 = tab(4);   % Oil viscosity at reference pressure.
   Cv  = tab(5);   % Oil viscosibility.

   dp = p - p0;

   B  = B0 .* exp(- C .* dp);
   dB = -C .* B;

   [R, dR] = deal(zeros(size(p)));

   y  = (C - Cv) .* dp;
   mu = mu0 .* B0 .* exp(-y) ./ B;  % == B*mu ./ B

   P = true(size(p));
end
