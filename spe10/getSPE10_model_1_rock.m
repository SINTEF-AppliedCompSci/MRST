function rock = getSPE10_model_1_rock()
%Define rock properties for Model 1 of tenth SPE CSP
%
% SYNOPSIS:
%   rock = getSPE10_model_1_rock
%
% PARAMETERS:
%   None.
%
% RETURNS:
%   rock - Rock structure having fields 'perm' and 'poro' for the entire
%          100-by-1-by-20 cell geometry of the first benchmark dataset of
%          the tenth SPE comparative solution project.  The permeability
%          distribution is a correlated geostatistically generated field
%          whereas the porosity is homogenous (phi=0.2) throughout the
%          formation.
%
% NOTE:
%   The permeability data is returned in strict SI units (metres squared).
%
% SEE ALSO:
%   `getSPE10rock`, `make_spe10_data`.

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

   data = load(fullfile(getDatasetPath('spe10'), 'model1_data'));
   rock = data.rock;
end
