function [bc, changed] = convertPressureBoundaryConditionsToFlux(G, state, bc)
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

    if isempty(bc)
        changed = false;
    else
        isDir = strcmpi(bc.type, 'pressure');
        changed = any(isDir);
        if any(isDir)
            % Convert Dirichlet boundary conditions to flux
            % boundary conditions for the transport
            dirFace = bc.face(isDir);
            q = sum(state.flux(dirFace, :), 2);
            sgn = 1 - 2*(G.faces.neighbors(dirFace, 2) == 0);
            bc.value(isDir) = sgn.*q;
            [bc.type{isDir}] = deal('resflux');
        end
    end
end
