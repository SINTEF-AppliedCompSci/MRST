function dp = nozzleValve(v, rho, D, dischargeCoeff, flowtype)
% Nozzle valve model

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

A = pi*(D/2).^2;
if nargin < 5
    flowtype = 'massRate';
end
switch flowtype
    case 'velocity'
        % do nothing
    case 'volumeRate'
        v = v./A;
    case 'massRate'
        v = v./(rho*A);
    otherwise
        error(['Unknown flow type: ', flowtype]);
end

dp = - (sign(value(v)).*rho.*v.^2)./(2*dischargeCoeff.^2);
end
