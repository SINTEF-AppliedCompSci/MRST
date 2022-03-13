function plotCoarseningStep(p, G, I1, I2, L, U, no, flag)                  %#ok<INUSD>
%Make a plot of the results of one step in a coarsening algorithm
%
% SYNOPSIS:
%   plotCoarseningStep(p, G, I1, I2, L, U, no)
%   plotCoarseningStep(p, G, I1, I2, L, U, no, flag)
%
% PARAMETERS:
%   G    - grid structure
%   p    - partition vector
%   I1   - indicator used to determine the size of the blocks. I1 should
%          be above a lower threshold for all blocks
%   I2   - indicator used to determine the flow through each block. I2
%          should be below an upper threshold for all blocks
%   L    - lower bound on indicator I1
%   U    - upper bound on indicator I2
%   no   - step number in the algorithm
%   flag - use indicator I2 rather than partition vector to color the
%          coarse blocks

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
targs = {'FontSize',10,'FontWeight','normal'};
clf
subplot(2,3,[1 4]),
if nargin<8
   plotCellData(G, mod(p, 17), 'EdgeColor', 'none');
else
   plotCellData(G, I2,'FaceAlpha',.5, 'EdgeColor', 'none');
   outlineCoarseGrid(G, p);
end
view(2), axis tight off;
title(sprintf('Step %d: %d blocks', no, max(p)),targs{:});

plural = {'', 's'};

subplot(2,3,[2 3]);
L = L * sum(I1.*G.cells.volumes) / G.cells.num;
bar(accumarray(p, I1.*G.cells.volumes));
hold on, plot([1 max(p)], [L L], 'r', 'LineWidth', 2); hold off
axis tight, set(gca,'ylim',[0 4*L]);


nviol = sum(accumarray(p, I1 .* G.cells.volumes) < L);
title(sprintf('Volume indicator: lower bound, %d violation%s', ...
              nviol, plural{1 + (nviol ~= 1)}),targs{:});

subplot(2,3,[5 6]);
U = U * sum(I2.*G.cells.volumes) / G.cells.num;
bar(accumarray(p, I2.*G.cells.volumes));
hold on, plot([1 max(p)], [U U], 'r', 'LineWidth', 2); hold off
axis tight, set(gca,'ylim',[0 2*U]);

nviol = sum(accumarray(p, I2 .* G.cells.volumes) > U);
title(sprintf('Flow indicator: upper bound, %d violation%s', ...
              nviol, plural{1 + (nviol ~= 1)}),targs{:});
