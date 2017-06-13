function state = updateAdsorption(state0, state, model)
%
%
% SYNOPSIS:
%   function state = updateAdsorption(state0, state, model)
%
% DESCRIPTION: Update the adsorption value in the state variable. Used by the
% surfactant models.
%
% PARAMETERS:
%   state0 - State at previous time step
%   state  - State at current time state
%   model  - Instance of the model
%
% RETURNS:
%   state - State at current time state with updated adsorption values.
%
% EXAMPLE:
%
% SEE ALSO: ExplicitConcentrationModel, FullyImplicitOilWaterSurfactantModel, ImplicitExplicitOilWaterSurfactantModel
%

%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

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


    c = model.getProps(state, 'surfactant');
    adsmax0 = model.getProps(state0, 'adsmax');

    ads = computeEffAds(c, adsmax0, model.fluid);
    adsmax = max(ads, adsmax0);
    state = model.setProp(state, 'ads', ads);
    state = model.setProp(state, 'adsmax', adsmax);

end


