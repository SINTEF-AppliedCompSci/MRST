function [p, solvetime] = reconstructPressure(partition, pressure, A, rhs)
% Solve reconstruction problem for multiscale methods

%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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
    D = formReconstructionMatrix(A, partition);
    
    tmp1 = warning('query','MATLAB:nearlySingularMatrix');
    tmp2 = warning('query','MATLAB:singularMatrix');
    warning('off','MATLAB:nearlySingularMatrix')
    warning('off','MATLAB:singularMatrix')
    t = tic();
    p = mldivide(D, rhs - (A - D)*pressure);
    solvetime = toc(t);
    warning(tmp1.state,'MATLAB:nearlySingularMatrix')
    warning(tmp2.state,'MATLAB:singularMatrix')
end
