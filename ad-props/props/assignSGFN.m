function f = assignSGFN(f, sgfn, reg)
    [f.krG, pcOG, f.krPts.g, hasPC] = getFunctions(f, sgfn, reg);
    if hasPC
        f.pcOG = pcOG;
    end
end

function [krG, pcOG, pts, hasPC] = getFunctions(f, SGFN, reg)
    [krG, pcOG] = deal(cell(1, reg.sat));
    
    pts = deal(zeros(reg.sat, 4));
    hasPC = false;
    for i = 1:reg.sat
        if isfield(f, 'krPts')
            swcon = f.krPts.w(i, 1);
        else
            warning('No relperm points found in assignment of SGOF.');
            swcon = 0;
        end
        pts(i, :) = getPoints(SGFN{i}, swcon);
        sgfn = extendTab(SGFN{i});
        SG = sgfn(:, 1);
        kr = sgfn(:, 2);
        pc = sgfn(:, 3);
        hasPC = hasPC || any(pc ~= 0);
        krG{i} = @(sg) interpTable(SG, kr, sg);
        pcOG{i} = @(sg) interpTable(SG, pc, sg);
    end
end

function pts = getPoints(sgfn, swcon)
    % Connate gas saturation
    pts = zeros(1, 4);
    pts(1) = sgfn(1, 1);
    % Last mobile gas saturation
    ii = find(sgfn(:,2)==0, 1, 'last');
    pts(2) = sgfn(ii,1);
    % Last point
    pts(3) = sgfn(end,1);
    % Maximum relperm
    pts(4) = sgfn(end,2);
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
