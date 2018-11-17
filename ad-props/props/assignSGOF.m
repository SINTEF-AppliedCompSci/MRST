function f = assignSGOF(f, sgof, reg)
    [f.krG, f.krOG, f.pcOG] = getFunctions(f, sgof, reg);
end

function [krG, krOG, pcOG] = getFunctions(f, SGOF, reg)
    [krG, krOG, pcOG] = deal(cell(1, reg.sat));
    
    for i = 1:reg.sat
        swof = extendTab(SGOF{i});
        SG = swof(:, 1);
        krG{i} = @(sg) interpTable(SG, swof(:, 2), sg);
        if isfield(f, 'sWcon')
            sWcon = f.sWcon(i);
        else
            sWcon = 0;
        end
        krOG{i} = @(so) interpTable(SG, swof(:, 3), 1-so-sWcon);
        pcOG{i} = @(sg) interpTable(SG, swof(:, 4), sg);
    end
end



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
