function f = assignPLYSHEAR(f, plyshear, reg)
% Polymer shear thinning/thickening
    f.plyshearMult = getFunction(plyshear, reg);
end

function plyshearMult = getFunction(plyshear, reg)
    plyshearMult = cell(1, reg.pvt);
    for i = 1:reg.pvt
        t = extendTab(plyshear{i});
        plyshearMult{i} = @(vW, varargin) reg.interp1d(t(:, 1), t(:, 2), vW);
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
