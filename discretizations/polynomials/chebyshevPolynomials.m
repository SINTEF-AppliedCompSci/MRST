function t = chebyshevPolynomials(degree)
%Undocumented Utility Function

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

    n = polyDim(degree, 1);
    t = cell(n,1);
    
    t{1} = Polynomial(0, 1);
    x    = Polynomial(1, 1);
    if degree > 0
        t{2} = Polynomial(1,1);
        for k = 1:n-2
            t{k+2} = 2*x*t{k+1} - t{k};
        end
    end
    
end
