function f = assignSGOF(f, sgof, reg)
    [f.krG, f.krOG, pcOG, f.krPts.g, f.krPts.og, hasPC] = getFunctions(f, sgof, reg);
    if hasPC
        f.pcOG = pcOG;
    end
end

function [krG, krOG, pcOG, pts, pts_o, hasPC] = getFunctions(f, SGOF, reg)
    [krG, krOG, pcOG] = deal(cell(1, reg.sat));
    
    [pts, pts_o] = deal(zeros(reg.sat, 4));
    hasPC = false;
    for i = 1:reg.sat
        if isfield(f, 'krPts')
            swcon = f.krPts.w(i, 1);
        else
            warning('No relperm points found in assignment of SGOF.');
            swcon = 0;
        end
        [pts(i, :), pts_o(i, :)] = getPoints(SGOF{i}, swcon);
        
        sgof = SGOF{i};
        SG = sgof(:, 1);
        ds = diff(SG);
        if isfield(reg, 'optimize') && reg.optimize && all(abs(ds - ds(1)) < sqrt(eps(ds(1))))
            % Uniform grid
            ds = ds(1);
            sgof = [[SG(1) - ds; SG; SG(end) + ds], sgof([1, 1:end, end], 2:end)];
            interp = reg.interp1d_uniform;
        else
            sgof = extendTab(sgof);
            interp = reg.interp1d;
        end
        SG = sgof(:, 1);
        pc = sgof(:, 4);
        hasPC = hasPC || any(pc ~= 0);
        SG_plus_swcon = SG(end:-1:1) + swcon;
        krog = sgof(end:-1:1, 3);
        krG{i} = @(sg) interp(SG, sgof(:, 2), sg);
        krOG{i} = @(so) interp(1-SG_plus_swcon, krog, so);
        pcOG{i} = @(sg) interp(SG, pc, sg);
    end
end

function [pts, pts_o] = getPoints(sgof, swcon)
    % Connate gas saturation
    pts = zeros(1, 4);
    pts(1) = sgof(1, 1);
    % Last mobile gas saturation
    ii = find(sgof(:,2)==0, 1, 'last');
    pts(2) = sgof(ii,1);
    % Last point
    pts(3) = sgof(end,1);
    % Maximum relperm
    pts(4) = sgof(end,2);
    
    % Get OG-scaling
    pts_o = zeros(1, 4);
    pts_o(3) = 1;
    ii = find(sgof(:,3) == 0, 1, 'first');
    pts_o(2) = 1 - sgof(ii,1) - swcon;
    % Maximum oil relperm
    pts_o(4) = sgof(1,3);
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
