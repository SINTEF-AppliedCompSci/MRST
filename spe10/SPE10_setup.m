function [G, W, rock] = SPE10_setup(varargin)
%Initialise properties for Model 2 of tenth SPE Comparative Solution Project
%
% This function is DEPRECATED and will be removed in a future version of
% MRST.  Please switch to using function getSPE10setup instead.
%
% SYNOPSIS:
%   [G, W, rock] = SPE10_setup
%   [G, W, rock] = SPE10_setup(layers)
%   [G, W, rock] = SPE10_setup(layers, wloc)
%
% PARAMETERS:
%   layers - Which of the 85 model layers to include in a specific test.
%            OPTIONAL.  If unspecified or empty, all 85 layers (a total of
%            60-by-220-by-85 == 1,122,000 cells) are included.
%
%            Some possible values are
%               layers = ( 1 : 35).';  %  Tarbert formation
%               layers = (36 : 85).';  %  Upper Ness formation
%
%   wloc   - Location of the five wells. OPTIONAL. If unspecified, we use
%            the default location
%               wloc     = [  1,   60,     1,   60,   30;
%                             1,    1,   220,  220,  110];
%
% RETURNS:
%   G    - MRST grid structure as described in grid_structure.
%
%   W    - Well structure.  Injector at 500 Bar, producers at 200 Bar.
%          Inner product 'ip_tpf'.
%
%   rock - Rock structure having fields 'perm' and 'poros' pertaining to
%          the specified layers.
%
% EXAMPLE:
%   [G, W, rock] = SPE10_setup(85)
%
% SEE ALSO:
%   `SPE10_rock`, `computeGeometry`.

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

   persistent WarningEmitted;

   if isempty(WarningEmitted),
      warning('SPEFunc:Deprecated', ...
             ['Function %s is deprecated and will be removed in a ', ...
              'future version of MRST\nPlease use the replacement ', ...
              'function getSPE10setup instead.'], mfilename);

      WarningEmitted = true;
   end

   [G, W, rock] = getSPE10setup(varargin{:});
end
