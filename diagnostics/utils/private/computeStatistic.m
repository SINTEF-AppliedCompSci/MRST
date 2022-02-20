function val = computeStatistic(vals, stat, prop)
%Undocumented Utility Function

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

type = 1; % standard mean/var
if nargin == 3
    if strcmp('TOF', prop(1:3))
        type = 2; % harmonic mean/var for TOFs
    end
end

n = size(vals,2);
switch stat
    case 'mean'
        if type == 1
            val = mean(vals, 2);
        else % harmonic
            val = 1./(mean(1./vals, 2));
        end
    case 'std'
        if type == 1
            val = std(vals, 0, 2);
        else % harmonic
            val = 1./(std(1./vals, 0, 2)); % or something ...
        end
    case 'max diff'
        val = max(vals, [], 2) - min(vals, [], 2);
end
end
