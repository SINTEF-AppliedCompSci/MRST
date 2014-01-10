function plotCoarseningStep(p, G, I1, I2, L, U, no, flag)
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
clf
subplot(2,3,[1 4]),
if nargin<8
   plotCellData(G, mod(p, 17));
else
   plotCellData(G, I2,'FaceAlpha',.5);
   h=outlineCoarseGrid(G, p);
end
view(2), axis tight off;
title(sprintf('Step %d: %d blocks', no, max(p)));

plural = {'', 's'};

subplot(2,3,[2 3]);
L = L * sum(I1.*G.cells.volumes) / G.cells.num;
bar(accumarray(p, I1.*G.cells.volumes));
hold on, plot([1 max(p)], [L L], 'r', 'LineWidth', 2); hold off
axis tight, set(gca,'ylim',[0 4*L]);


nviol = sum(accumarray(p, I1 .* G.cells.volumes) < L);
title(sprintf('Volume indicator: lower bound, %d violation%s', ...
              nviol, plural{1 + (nviol ~= 1)}));

subplot(2,3,[5 6]);
U = U * sum(I2.*G.cells.volumes) / G.cells.num;
bar(accumarray(p, I2.*G.cells.volumes));
hold on, plot([1 max(p)], [U U], 'r', 'LineWidth', 2); hold off
axis tight, set(gca,'ylim',[0 2*U]);

nviol = sum(accumarray(p, I2 .* G.cells.volumes) > U);
title(sprintf('Flow indicator: upper bound, %d violation%s', ...
              nviol, plural{1 + (nviol ~= 1)}));
