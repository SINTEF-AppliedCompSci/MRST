function kr_deck = getSPE10_model_1_relperm()
%Define Oil/Gas Relative Permeability Properties for Model 1 of tenth SPE CSP
%
% SYNOPSIS:
%   kr_deck = getSPE10_model_1_relperm
%
% PARAMETERS:
%   None.
%
% RETURNS:
%   kr_deck - Stripped-down simulation deck for ADI-fluid (relperm)
%             construction.  This input-deck contains just enough basic
%             information to form relative permeability curves from the CSP
%             benchmark data through the use of module 'ad-props's function
%             initDeckADIFluid.
%
% SEE ALSO:
%   `getSPE10_model_1_rock`, `initDeckADIFluid`.

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

   data    = load(fullfile(getDatasetPath('spe10'), 'model1_data'));
   kr_deck = data.kr_deck;
end
