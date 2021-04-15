function f = assignSWFN(f, swfn, reg)
    [f.krW, pcOW, f.krPts.w, hasPC] = getFunctions(swfn, reg);
    if hasPC
        f.pcOW = pcOW;
    end

    if isfield(reg, 'SURFNUM')
       % Assign miscible relperm for surfactant
       f.krWSft  = @(sw, varargin)krWSft(sw, swfn, reg, varargin{:});
       % Assign residual water saturation for surfactant
       swcon = f.krPts.w(:, 1);
       f.sWconSft = swcon(reg.SURFNUM);
       % Assign residual oil saturation
       sOres  = cellfun(@(x)x(end, 1), swfn);
       f.sOres = 1 - sOres(reg.SATNUM);
       f.sOresSft = 1 - sOres(reg.SURFNUM);
    end

end


function [krW, pcOW, pts_w, hasPC] = getFunctions(SWFN, reg)
    pts_w = zeros(reg.sat, 4);
    
    [krW, pcOW] = deal(cell(1, reg.sat));
    hasPC = false;
    for i = 1:reg.sat
        pts_w(i, :) = getPoints(SWFN{i});
        
        swfn = extendTab(SWFN{i});
        SW = swfn(:, 1);
        kr = swfn(:, 2);
        pc = swfn(:, 3);
        hasPC = hasPC || any(pc ~= 0);
        krW{i} = @(sw) interpTable(SW, kr, sw);
        pcOW{i} = @(sw) interpTable(SW, pc, sw);
    end
end


function pts = getPoints(swof)
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
end

function v = krWSft(sw, swfn, reg, varargin)
surfinx = getRegMap(sw, reg.SURFNUM, reg.SURFINX, varargin{:});
T = cellfun(@(x)x(:,[1,2]), swfn, 'UniformOutput', false);
T = extendTab(T);
v = interpReg(T, sw, surfinx);
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
