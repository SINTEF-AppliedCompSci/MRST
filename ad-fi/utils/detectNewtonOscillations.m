function [isOscillating, isStagnating] = detectNewtonOscillations(history, primary, current, tol)
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

    if current < 3
        isOscillating = false;
        isStagnating = false;
        return
    end

    tmp = history(current-2:current, primary);

    oscillate =  relChange(tmp(1,:), tmp(3,:)) < tol & ...
                 relChange(tmp(2,:), tmp(3,:)) > tol;

    stagnate = relChange(tmp(2,:), tmp(1,:));
    stagnate(isnan(stagnate)) = 0;

    isStagnating = all(stagnate < 1e-3);
    isOscillating = sum(oscillate) > 1;
end

function v = relChange(a,b)
    v = abs((a-b)./b);
end
