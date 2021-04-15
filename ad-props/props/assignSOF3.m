function f = assignSOF3(f, sof3, reg)
    [f.krOW, f.krOG, f.krPts.ow, f.krPts.og] = getFunctions(sof3, reg);

end

function [krOW, krOG, pts_ow, pts_og] = getFunctions(SOF3, reg)
    [krOW, krOG] = deal(cell(1, reg.sat));
    
    [pts_ow, pts_og] = deal(zeros(reg.sat, 4));
    for i = 1:reg.sat
        sof3 = SOF3{i};
        SO = sof3(:, 1);

        pts_ow(i, :) = getPoints(SO, sof3(:, 2));
        pts_og(i, :) = getPoints(SO, sof3(:, 3));
        sof3 = extendTab(sof3);
        SO = sof3(:, 1);
        krow = sof3(:, 2);
        krog = sof3(:, 3);
        krOW{i} = @(so) interpTable(SO, krow, so);
        krOG{i} = @(so) interpTable(SO, krog, so);
    end
end

function pts = getPoints(so, kro)
    pts = zeros(1, 4);
    ii = find(kro == 0, 1, 'first');
    pts(2) = so(ii);
    pts(3) = 1;
    pts(4) = kro(end);
end



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

