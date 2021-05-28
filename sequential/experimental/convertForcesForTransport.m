function forces = convertForcesForTransport(state, forces)
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

    if ~isfield(forces, 'bc') || isempty(forces.bc)
        return
    end
    isDir = strcmpi(forces.bc.type, 'pressure');
    [forces.bc.type{isDir}] = deal('flux');
    
    qBC = state.flux(forces.bc.face(isDir), :);
    qT = sum(qBC, 2);
    netq = sum(abs(qBC), 2);
    
    sat0 = forces.bc.sat;
    forces.bc.value(isDir) = qT;
    forces.bc.sat(isDir, :) = bsxfun(@rdivide, abs(qBC), netq);

    bad = any(isnan(forces.bc.sat), 2);
    forces.bc.sat(bad, :) = sat0(bad, :);
end
