function rock = SPE10_rock(varargin)
%Define rock properties for Model 2 of tenth SPE CSP
%
% This function is DEPRECATED and will be removed in a future version of
% MRST.  Please switch to using function getSPE10rock instead.
%
% SYNOPSIS:
%   rock = SPE10_rock
%   rock = SPE10_rock(layers)
%   rock = SPE10_rock(I, J, K)
%
% PARAMETERS:
%   layers - Which of the 85 model layers to include in a specific test.
%            OPTIONAL.  If unspecified, all 85 layers (for a total of
%            60-by-220-by-85 == 1,122,000 cells) are included.
%
%            Some possible values are
%               layers = ( 1 : 35).';  %  Tarbert formation
%               layers = (36 : 85).';  %  Upper Ness formation
%
%   I,J,K  - Global Cartesian indices identifying subset of rock data.
%            Useful, e.g., in extracting lateral sections.  Arrays must
%            satisfy
%
%               ALL((0 < I) & (I <=  60)) && ...
%               ALL((0 < J) & (J <= 220)) && ...
%               ALL((0 < K) & (K <=  85))
%
%            OPTIONAL.  If unspecified, treated as I=1:60,J=1:220,K=1:85
%            (i.e., the entire model; 1,122,000 cells).  In other words,
%            equivalently to an unspecified 'layers' input.
%
% RETURNS:
%   rock - Rock structure having fields 'perm' and 'poro' pertaining to
%          the specified layers or cell subset.
%
% NOTE:
%   The permeability data is returned in units of milli*darcy.
%
%   It is the caller's responsibility to convert this data into MRST's
%   strict SI-only unit conventions.  Function 'convertFrom' provides one
%   conversion method.
%
% EXAMPLE:
%   rock = SPE10_rock(85)
%
% SEE ALSO:
%   `getSPE10rock`, `SPE10_setup`, `convertFrom`, `milli`, `darcy`.

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
             ['Function %s is deprecated and will be removed in a ' , ...
              'future version of MRST.\nPlease use the replacement ', ...
              'function getSPE10rock instead.\n'                    , ...
              'Please make special note of changed unit convention ', ...
              'for permeabilities in getSPE10rock.'], mfilename);

      WarningEmitted = true;
   end

   rock      = getSPE10rock(varargin{:});
   rock.perm = convertTo(rock.perm, milli*darcy);
end
