function h = colorizeCatchmentRegions( Gt, ta )
% Add catchment areas to current figure window

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

    % Plot accumulation regions
    if max(ta.traps)<128
       cmap = lines(max(ta.traps) + 1);
    else
       cmap = jet(max(ta.traps) + 1);
    end
    colormap(cmap);
    %map = greyMask(cmap);
    map = (cmap + repmat(get(gcf, 'Color'), size(cmap, 1), 1))./2;
    map(1,:) = get(gcf, 'Color');
    tmp = ta.trap_regions;
    tmp(tmp>max(ta.traps)) = max(ta.traps);

    h = plotCellData(Gt, ones(Gt.cells.num, 1), 'EdgeColor', 'none');
    set(h, 'FaceVertexCData', map(tmp + 1, :))
end
