function f = assignGSF(f, gsf, reg)
%Process Saturation Functions for Gas in Gas/Water System (GSF Keyword)
%
% SYNOPSIS:
%   f = assignGSF(f, wsf, reg)
%
% PARAMETERS:
%   f   - Fluid structure.
%
%   gsf - Cell array of m-by-r GSF tables.  The first column is the
%         gas saturation, the second column is the relative
%         permeability of gas, and the third column is the gas/water
%         capillary pressure (i.e., Pg - Pw).
%
%   reg - Saturation function region description.
%
% RETURNS:
%   f - Updated fluid structure, now including saturation functions
%       for the gas phase.

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

   [f.krG, pcWG, f.krPts.g, hasPC] = getFunctions(gsf, reg);

   if hasPC
      f.pcWG = pcWG;
   end
end

% --------------------------------------------------------------------------

function [krG, pcWG, pts, hasPC] = getFunctions(GSF, reg)
   [krG, pcWG] = deal(cell([1, reg.sat]));
   pts   = zeros([reg.sat, 4]);
   hasPC = false;

   for r = 1 : reg.sat
      pts(r, :) = getPoints(GSF{r});

      [krG{r}, pcWG{r}, nonZeroPc] = ...
         regionInterpolants(extendTab(GSF{r}));

      if nonZeroPc
         hasPC = true;
      end
   end
end

% --------------------------------------------------------------------------

function [krG, pcWG, nonZeroPc] = regionInterpolants(gsf);
   krG  = @(sg) interpTable(gsf(:,1), gsf(:,2), sg);
   pcWG = @(sg) interpTable(gsf(:,1), gsf(:,3), sg);

   nonZeroPc = any(gsf(:,3) ~= 0);
end

% --------------------------------------------------------------------------

function pts = getPoints(gsf)
   pts = zeros([1, 4]);

   % Connate gas saturation.
   pts(1) = gsf(1, 1);

   % Last immobile gas saturation.
   ix = find(gsf(:,2) > 0, 1, 'first');
   pts(2) = gsf(ix - 1, 1);

   % Maximum gas saturation.
   pts(3) = gsf(end, 1);

   % Maximum relperm.
   pts(4) = gsf(end, 2);
end
