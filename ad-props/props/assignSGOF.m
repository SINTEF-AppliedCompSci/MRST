function f = assignSGOF(f, sgof, reg)
    [f.krG, f.krOG, f.pcOG, f.krPts.g, f.krPts.og] = getFunctions(f, sgof, reg);
end

function [krG, krOG, pcOG, pts, pts_o] = getFunctions(f, SGOF, reg)
    [krG, krOG, pcOG] = deal(cell(1, reg.sat));
    
    [pts, pts_o] = deal(zeros(reg.sat, 3));
    for i = 1:reg.sat
        [pts(i, :), pts_o(i, :)] = getPoints(SGOF{i});
        swof = extendTab(SGOF{i});
        SG = swof(:, 1);
        krG{i} = @(sg) interpTable(SG, swof(:, 2), sg);

        krOG{i} = @(so) interpTable(SG, swof(:, 3), 1-so-f.krPts.w(i, 1));
        pcOG{i} = @(sg) interpTable(SG, swof(:, 4), sg);
    end
end

function [pts, pts_o] = getPoints(sgof)
    % Connate gas saturation
    pts = zeros(1, 3);
    pts(1) = sgof(1, 1);
    % Last mobile gas saturation
    ii = find(sgof(:,2)==0, 1, 'last');
    pts(2) = sgof(ii,1);
    % Last point
    pts(3) = sgof(end,1);
    
    % Get OG-scaling
    pts_o = zeros(1, 3);
    pts_o(3) = 1;
    ii = find(sgof(:,3) == 0, 1, 'first');
    pts_o(2) = 1 - sgof(ii,1);
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
