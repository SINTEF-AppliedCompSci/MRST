function f = assignSWOF(f, swof, reg)
    [f.krW, f.krOW, pcOW, f.krPts.w, f.krPts.ow, hasPC] = getFunctions(swof, reg);
    if hasPC
        f.pcOW = pcOW;
    end
end

function [krW, krOW, pcOW, pts_w, pts_ow, hasPC] = getFunctions(SWOF, reg)
    pts_w = zeros(reg.sat, 4);
    pts_ow = zeros(reg.sat, 4);
    
    [krW, krOW, pcOW] = deal(cell(1, reg.sat));
    hasPC = false;
    for i = 1:reg.sat
        [pts_w(i, :), pts_ow(i, :)] = getPoints(SWOF{i});
        swof = SWOF{i};
        SW = swof(:, 1);
        ds = diff(SW);
        if reg.optimize && all(abs(ds - ds(1)) < sqrt(eps(ds(1))))
            % Uniform grid
            ds = ds(1);
            swof = [[SW(1) - ds; SW; SW(end) + ds], swof([1, 1:end, end], 2:end)];
            interp = reg.interp1d_uniform;
        else
            swof = extendTab(swof);
            interp = reg.interp1d;
        end

        SW = swof(:, 1);
        pc = swof(:, 4);
        hasPC = hasPC || any(pc ~= 0);
        % Reverse krO table to be increasing and given with respect to so
        krow = swof(end:-1:1, 3);
        SO = 1 - SW(end:-1:1);
        krW{i} = @(sw) interp(SW, swof(:, 2), sw);
        krOW{i} = @(so) interp(SO, krow, so);
        pcOW{i} = @(sw) interp(SW, pc, sw);
    end
end


function [pts, pts_o] = getPoints(swof)
    pts = zeros(1, 4);
    % Connate water saturation
    pts(1) = swof(1, 1);
    % Last immobile water saturation
    ii = find(swof(:,2)==0, 1, 'last');
    pts(2) = swof(ii,1);
    % Last point
    pts(3) = swof(end,1);
    % Maximum relperm
    pts(4) = swof(end,2);
    
    % Get OW-scaling
    pts_o = zeros(1, 4);
    pts_o(3) = 1;
    ii = find(swof(:,3) == 0, 1, 'first');
    pts_o(2) = 1 - swof(ii,1);
    % Maximum oil relperm
    pts_o(4) = swof(1,3);
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

