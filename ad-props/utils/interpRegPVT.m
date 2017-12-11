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

if nreg > 1
    reginxSat  = cellfun(@(v) v(flag(v)), reginx, 'UniformOutput', false);
    reginxUSat = cellfun(@(v) v(~flag(v)), reginx, 'UniformOutput', false);
else
    reginxSat{1}  = ':';
    reginxUSat{1} = ':';
end

if ~compDer
    yi(flag)  = interpReg(tabSat, xi(flag), reginxSat);
    yi(~flag) = interp2DPVT(T, xi(~flag), vi(~flag), reginxUSat);
else
    [yi(flag), dyidxi(flag)] = interpReg(tabSat, xi(flag), reginxSat);
    [yi(~flag), dyidxi(~flag), dyidvi(~flag)] = interp2DPVT(T, xi(~flag), vi(~flag), reginxUSat);
end

end


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

