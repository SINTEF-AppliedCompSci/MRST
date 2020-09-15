function [yi, dyidxi, dyidvi] = interpPVT(T, xi, vi, flag)
% Interpolate PVT-type curves
compDer = (nargout>1);
if numelValue(flag) == 1
    flag = repmat(flag, numelValue(xi), 1);
end
if isa(T, 'function_handle')
    if compDer
        [yi, dyidxi, dyidvi] = T(xi, vi, flag);
    else
        yi = T(xi, vi, flag);
    end
else
    yi = zeros(size(xi));
    if compDer
        dyidxi = zeros(size(xi));
        dyidvi = zeros(size(xi));
    end
    tabSat = T.data(T.pos(1:end-1),:);

    if isempty(xi)
        [yi, dyidxi, dyidvi] = deal([]);
        return;
    else
        [ixf, ixnf] = deal(flag, ~flag);
        if ~compDer
            yi(ixf)  = interpTable(tabSat(:, 1), tabSat(:, 2), xi(ixf));
            yi(ixnf) = interp2DPVT({T}, xi(ixnf), vi(ixnf), {':'});
        else
            [yi(ixf), dyidxi(ixf)]                 = interpReg({tabSat}, xi(ixf), {':'});
            [yi(ixnf), dyidxi(ixnf), dyidvi(ixnf)] = interp2DPVT({T}, xi(ixnf), vi(ixnf), {':'});
        end
    end
end
%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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

