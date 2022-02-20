function varargout = plotWellPairConnections(G, WP, D, W, pv, minAlloc)
%Plot lines between wells to show relative flux allocation
%
% SYNOPSIS
%   plotWellPairConnections(G, WP, D, W, pv)
%   plotWellPairConnections(G, WP, D, W, pv, mA)
%
% PARAMETERS:
%   G     - Grid structure.
%
%   WP - data structure containing information about well pairs,
%        computed by a call to 'computeWellPairs'
%
%   D  - data structure with basic data for flow diagnostics computed
%        by a call to 'computeTOFandTracer'
%
%   W  - Well structure as defined by function 'addWell'.
%
%   pv - vector with pore volume as computed by 'poreVolume(G,rock)'
%
%   mA - disregard connections with relative allocation below mA
%
% DESCRIPTION:
%   The routine makes a curved line between each well pair, where the
%   thickness of the line is proportional to the relative flux allocation
%   for this pair.
%
% SEE ALSO:
%   `computeTOFandTracer`, `computeWellPairs`
%
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

if nargin<6
    minAlloc = 0.01;
end
hold on

maxAlloc = 0;
for i = 1:numel(WP.inj)
    maxAlloc = max([maxAlloc sum(WP.inj(i).alloc, 1)]);
end

handles = [];
for i=1:numel(D.inj)
    ipos = mean(G.cells.centroids(W(D.inj(i)).cells,:), 1);
    ialloc = sum(WP.inj(i).alloc, 1);

    for p=1:numel(D.prod)
        ppos = mean(G.cells.centroids(W(D.prod(p)).cells,:), 1);
        alloc = ialloc(p);

        if alloc/sum(ialloc) > minAlloc
            localRegion = find(D.ipart == i & D.ppart == p);
            center = sum(G.cells.centroids(localRegion, :).* ...
                repmat(pv(localRegion), 1, G.griddim), 1)/sum(pv(localRegion));

            pts = [ipos; center; ppos];

            thick = 20*alloc/(maxAlloc);

            htext = text(center(1), center(2), center(3)-2*thick, ...
                sprintf('%2.1f%%', 100*(alloc/sum(ialloc))),...
                'VerticalAlignment', 'Bottom', 'FontWeight', ...
                'bold', 'FontSize', 12, 'Color', 'k');

            hline = plot3(pts(:,1), pts(:,2), pts(:,3),...
                '-','LineWidth', thick, 'Color', colorizeWell('prod', p, D));
            handles = [handles; hline; htext];                             %#ok<AGROW>
        end
    end
end

if nargout > 0
    varargout{1} = handles;
end
hold off
end
