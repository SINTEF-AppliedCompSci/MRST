function st = getCellStatusVO(model, state, sO, sW, sG)
% Get status flags for each cell in a black-oil model
%
% SYNOPSIS:
%   st = getCellStatusVO(model, state, sO, sW, sG)
%
% DESCRIPTION:
%   Get the status flags for the number of phases present.
%
% REQUIRED PARAMETERS:
%   model     - ThreePhaseBlackOilModel-derived class. Typically,
%               equationsBlackOil will be called from the class
%               getEquations member function.
%
%   state     - The state for which the status flags are to be computed.
%
%   sO        - Oil saturation. One value per cell in the simulation grid.
%
%   sW        - Water saturation. One value per cell in the simulation grid.
%
%   sG        - Gas saturation. One value per cell in the simulation grid.
%
% RETURNS:
%   st - Cell array with three columns with one entry per cell. The
%        interpretation of each flag: 
%      Col 1: A cell is flagged as true if oil is present, but gas
%      is not present.
%      Col 2: A cell is flagged as true if gas is present, but oil
%      is not present.
%      Col 3: A cell is flagged as true if both gas and oil are present and
%      true three-phase flow is occuring.
%
% SEE ALSO:
%   ThreePhaseBlackOilModel

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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

    if isfield(state, 'status')
        % Status should be passed on from updateStateVO (to be sure definition is
        % identical). rs and rv are assumed to be compatible, i.e. rx = rxSat for
        % saturated cells and rx <= rxSat for undersaturated. Three values of
        % status are:
        % status 0: should not occur (almost water only -> state 3)
        % status 1 oil, no gas  : x = rs, sg = 0    , rv = rvMax
        % status 2 gas, no oil  : x = rv, sg = 1-sw , rs = rsMax
        % status 3 oil and gas  : x = sg, rs = rsMax, rv = rvMax
        status = state.status;
    else
        watOnly    = sW > 1- sqrt(eps);
        if ~model.vapoil
            oilPresent = true;
        else
            oilPresent = or(sO > 0, watOnly);
        end
        if ~model.disgas
            gasPresent = true;
        else
            gasPresent = or(sG > 0, watOnly);
        end
        status = oilPresent + 2*gasPresent;
    end

    if ~model.disgas
        st1 = false;
    else
        st1 = status==1;
    end
    if ~model.vapoil
        st2 = false;
    else
        st2 = status==2;
    end
    st3 = status == 3;
    st = {st1, st2, st3};
end