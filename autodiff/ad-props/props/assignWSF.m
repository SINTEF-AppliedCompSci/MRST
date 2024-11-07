function f = assignWSF(f, wsf, reg)
%Process Saturation Functions for Water in Gas/Water System (WSF Keyword)
%
% SYNOPSIS:
%   f = assignWSF(f, wsf, reg)
%
% PARAMETERS:
%   f   - Fluid structure.
%
%   wsf - Cell array of m-by-2 WSF tables.  The first column is the
%         water saturation while the second column is the relative
%         permeability of water.
%
%   reg - Saturation function region description.
%
% RETURNS:
%   f - Updated fluid structure, now including saturation functions
%       for the water phase.

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

   [f.krW, f.krPts.w] = getFunctions(wsf, reg);
end

% --------------------------------------------------------------------------

function [krW, pts] = getFunctions(WSF, reg)
   pts = zeros([reg.sat, 4]);
   krW = cell([1, reg.sat]);

   for r = 1 : reg.sat
      pts(r, :) = getPoints(WSF{r});
      krW{r}    = waterRelpermInterpolant(WSF{r}, reg);
   end
end

% --------------------------------------------------------------------------

function krW = waterRelpermInterpolant(wsf, reg)
   SW = wsf(:, 1);
   ds = diff(SW);

   if reg.optimize && all(abs(ds - ds(1)) < sqrt(eps(ds(1))))
      % Uniform grid
      ds = ds(1);
      wsf = [[SW(1) - ds; SW; SW(end) + ds], wsf([1, 1:end, end], 2)];
      interp = reg.interp1d_uniform;
   else
      wsf = extendTab(wsf);
      interp = reg.interp1d;
   end

   krW = @(sw) interp(wsf(:,1), wsf(:,2), sw);
end

% --------------------------------------------------------------------------

function pts = getPoints(wsf)
   pts = zeros([1, 4]);

   % Connate water saturation
   pts(1) = wsf(1, 1);

   % Last immobile water saturation
   ii = find(wsf(:,2) > 0, 1, 'first');
   pts(2) = wsf(ii - 1, 1);

   % Last point
   pts(3) = wsf(end, 1);

   % Maximum relperm
   pts(4) = wsf(end, 2);
end
