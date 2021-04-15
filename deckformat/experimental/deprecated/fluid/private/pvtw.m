function [B, dB, R, dR, mu, P] = pvtw(tab, p)
%Evaluate water PVT properties.
%
% SYNOPSIS:
%   [B, dB, R, dR, mu, P] = pvtw(tab, p)
%
% PARAMETERS:
%
% RETURNS:
%   B, dB - Volume formation factor (B) and derivative of B with respect to
%           pressure (dB).
%
%   R, dR - Miscibility ratio (R) and derivative of R with respect to
%           pressure (dR).  Water does not mix, so both R and dR == 0.
%
%   mu    - Water viscosity.
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
   B0  = tab(2);   % Water formation factor at p0.
   Cw  = tab(3);   % Water compressibility.
   mu0 = tab(4);   % Water viscosity at p0.
   Cv  = tab(5);   % Water viscosibility.

   if abs(Cw) > 0,
      if(true)
         x  = Cw .* (p - p0);
         d  = 1 + x + x.^2./2;
         B  =   B0 ./ d;
         dB = - B0 ./ d.^2 .* (1 + x) .* Cw;
      else
         B = B0.*exp(- Cw.*(p-p0));
         dB = -Cw.*B;
      end
   else
      % Incompressible.
      B  = ones(size(p));
      dB = zeros(size(p));
   end

   [R, dR] = deal(zeros(size(p)));

   if abs(Cv) > 0,
      % E300 formula.
      y  = -Cv .* (p - p0);
      mu = mu0 ./ (1 + y + y.^2./2);
   else
      % Incompressible.
      mu = repmat(mu0, size(p));
   end

   P = true(size(p));
end
