function [sG, rs, rv, rsSat, rvSat] = calculateHydrocarbonsFromStatusBO(fluid, ...
                                                      status, sO, x, rs, rv, ...
                                                      pressure, disgas, vapoil)
% Compute solution variables for the gas/oil/rs/rv-variable in black-oil
%
% SYNOPSIS:
%   [sG, rs, rv, rsSat, rvSat] = ...
% calculateHydrocarbonsFromStatusBO(model, status, sO, x, rs, rv, pressure)
%
% DESCRIPTION:
%   The ThreePhaseBlackOil model has a single unknown that represents
%   either gas, oil, oil in gas (rs) and gas in oil (rv) on a cell-by-cell
%   basis. The purpose of this function is to easily compute the different
%   quantities from the "x"-variable and other properties, with correct
%   derivatives if any of them are AD-variables.
%
% REQUIRED PARAMETERS:
%   model    - ThreePhaseBlackOilModel-derived class. Determines if
%              vapoil/disgas is being used and contains the fluid model.
%
%   status   - Status flag as defined by "getCellStatusVO"
%
%   sO       - The tentative oil saturation.
%
%   x        - Variable that is to be decomposed into sG, sO, rs, rv, ...
%
%   rs, rv   - Dissolved gas, vaporized oil
%
%   pressure - Reservoir oil pressure
%
%   disgas   - true if dissolved gas should be taken into account
%
%   vapoil   - true if vaporized oil should be taken into account
%
% OPTIONAL PARAMETERS:
%   'field'   -  
% RETURNS:
%   
%
% SEE ALSO:
%   equationsBlackOil, getCellStatusVO

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



    sG = status{2}.*sO + status{3}.*x;
    if disgas
        if iscell(fluid.rsSat)
            rsSat = fluid.rsSat{1}(pressure);
        else
            rsSat = fluid.rsSat(pressure);
        end
        rs = (~status{1}).*rsSat + status{1}.*x;
    else % otherwise rs = rsSat = const
        rsSat = rs;
    end
    if vapoil
        % Note: This is a hack.
        if iscell(fluid.rvSat)
            rvSat = fluid.rvSat{1}(pressure);
        else
            rvSat = fluid.rvSat(pressure);
        end
        rv = (~status{2}).*rvSat + status{2}.*x;
    else % otherwise rv = rvSat = const
        rvSat = rv;
    end
end

