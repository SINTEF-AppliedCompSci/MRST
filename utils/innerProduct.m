function ip = innerProduct(u, v, cubature, elements, type)
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

    if isa(u, 'SpatialVector') || isa(v, 'SpatialVector')
        uv = dot(u, v);
    else
        uv = u.*v;
    end
    switch type
        case 'cell'
            ip = cellIntegral(uv, cubature, elements);
        case 'face'
            ip = faceIntegral(uv, cubature, elements);
    end
end
