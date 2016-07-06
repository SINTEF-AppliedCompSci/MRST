function [sG, rs, rv, rsSat, rvSat] = calculateHydrocarbonsFromStatusBO(model, status, sO, x, rs, rv, pressure)
    % define sG, rs and rv in terms of x
    % We cast "sO" to double to avoid cancellation of derivatives when
    % calculating the final sO later.
    sG = status{2}.*double(sO) + status{3}.*x;
    if model.disgas
        rsSat = model.fluid.rsSat(pressure);
        rs = (~status{1}).*rsSat + status{1}.*x;
    else % otherwise rs = rsSat = const
        rsSat = rs;
    end
    if model.vapoil
        rvSat = model.fluid.rvSat(pressure);
        rv = (~status{2}).*rvSat + status{2}.*x;
    else % otherwise rv = rvSat = const
        rvSat = rv;
    end
end

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
