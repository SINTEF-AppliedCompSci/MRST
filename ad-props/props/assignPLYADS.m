function f = assignPLYADS(f, plyads, reg)
   plyads = extendTab(plyads);
   f.ads = getFunction(plyads, reg);
end

function plyads = getFunction(PLYADS, reg)
   plyads = cell(1, reg.sat);
    for i = 1:reg.sat
        t = PLYADS{i};
        t = extendTab(t);
        plyads{i} = @(c, varargin) reg.interp1d(t(:, 1), t(:, 2), c);
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
