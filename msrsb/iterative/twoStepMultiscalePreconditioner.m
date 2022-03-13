function x = twoStepMultiscalePreconditioner(A, b, solveCoarse, smooth, it)
% Apply two step multiscale preconditioner. Internal routine.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
    if nargin < 5
        it = 1;
    end
    if it == 1
        x = smooth(b);
        x = x + solveCoarse(b - A*x);
    else
        x = zeros(size(b));
        for i = 1:it
            x = x + smooth(b - A*x);
            x = x + solveCoarse(b - A*x);
        end
    end
end
