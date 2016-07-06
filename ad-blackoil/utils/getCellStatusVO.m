function st = getCellStatusVO(model, state, oil, wat, gas)
% Status should be passed on from updateStateVO (to be sure definition is
% identical). rs and rv are assumed to be compatible, i.e. rx = rxSat for
% saturated cells and rx <= rxSat for undersaturated. Three values of
% status are:
% status 0: should not occur (almost water only -> state 3)
% status 1 oil, no gas  : x = rs, sg = 0    , rv = rvMax
% status 2 gas, no oil  : x = rv, sg = 1-sw , rs = rsMax
% status 3 oil and gas  : x = sg, rs = rsMax, rv = rvMax

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
        status = state.status;
    else
        watOnly    = wat > 1- sqrt(eps);
        if ~model.vapoil
            oilPresent = true;
        else
            oilPresent = or(oil > 0, watOnly);
        end
        if ~model.disgas
            gasPresent = true;
        else
            gasPresent = or(gas > 0, watOnly);
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