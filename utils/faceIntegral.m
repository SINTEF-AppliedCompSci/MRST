function int_u = faceIntegral(u, cubature, faces)
    % Integrate integrand u over cells using a given cubature
    % int_u(i) = (int_{face(i)} u ds)/|face(i)|

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

    if nargin < 3 || isinf(faces)
        % Empty cells means all internal connections of the grid
        faces = find(disc.internalConn)';
    end
    % Get cubature for all cells, transform coordinates to ref space
    W = cubature.getCubature(faces, 'face');
    int_u = W*u;
end
