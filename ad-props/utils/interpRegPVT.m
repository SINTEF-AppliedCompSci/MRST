function [yi, dyidxi, dyidvi] = interpRegPVT(T, xi, vi, flag, reginx)
% Interpolate PVT-type curves with region support
compDer = (nargout>1);
nreg = numel(reginx);

yi = zeros(size(xi));
if compDer
    dyidxi = zeros(size(xi));
    dyidvi = zeros(size(xi));
end

tabSat = cellfun(@(x)x.data(x.pos(1:end-1),:), T, 'UniformOutput', false);

if isempty(xi)
    [yi, dyidxi, dyidvi] = deal([]);
    return;
else
    for k = 1:nreg
        if nreg == 1
            [ixf, ixnf] = deal(flag, ~flag);
        else
            [ixf, ixnf] = deal(flag & reginx{k}, (~flag) & reginx{k});
        end
        if ~compDer
            yi(ixf)  = interpReg({tabSat{k}}, xi(ixf), {':'});
            yi(ixnf) = interp2DPVT({T{k}}, xi(ixnf), vi(ixnf), {':'});
        else
            [yi(ixf), dyidxi(ixf)]                 = interpReg({tabSat{k}}, xi(ixf), {':'});
            [yi(ixnf), dyidxi(ixnf), dyidvi(ixnf)] = interp2DPVT({T{k}}, xi(ixnf), vi(ixnf), {':'});
        end
    end
end


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
