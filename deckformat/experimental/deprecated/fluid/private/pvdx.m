function [B, dB, R, dR, mu, P] = pvdx(tab, p, varargin)
%Interpolate dead oil/dead gas PVT properties.
%
% SYNOPSIS:
%   [B, dB, R, dR, mu, P] = pvdx(tab, p)
%
% PARAMETERS:
%
% RETURNS:
%   B, dB - Volume formation factor (B) and derivative of B with respect to
%           pressure (dB).
%
%   R, dR - Miscibility ratio (R) and derivative of R with respect to
%           pressure (dR) (== 0 for oil/gas with no miscibility).
%
%   mu    - Phase viscosity.
%
%   P     - Predicate for identifying which cells are in a saturated state.
%
% SEE ALSO:
%   `pvto`, `pvtg`.

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


   % Interpolate 1/B
   [mu, B, dB] = interpViscVol(tab(:,1), tab(:, [3, 2]), p);
   R  = zeros([numel(p),1]);
   dR = zeros([numel(p),1]);
   P  = true(size(p));
end

%--------------------------------------------------------------------------

function [mu, B, dB] = interpViscVol(p, t, pi)
   if isempty(pi),
      [mu, B, dB] = deal([]);
   elseif numel(p) == 1,
      % Incompressible
      mu = repmat(t(1), [numel(pi), 1]);
      B  = repmat(t(2), [numel(pi), 1]);
      dB = zeros([numel(pi), 1]);
   else
      mu =            interpTable(p,      t(:,1), pi);
      B  =   1    ./  interpTable(p, 1 ./ t(:,2), pi);
      dB = - B.^2 .* dinterpTable(p, 1 ./ t(:,2), pi);
   end
end
