function myPlotPartition(G, W, p0, n)
% Plots partition, wells, and distribution of block volumes
%
% SYNOPSIS:
%     myPlotPartition(G, W, p, n)
%
% INPUT PARAMETERS:
%     G  - grid structure
%     W  - well structure
%     p  - partition vector
%     n  - color axis is set to [1,n] (optional)

persistent rp;
if nargin<4,
    n = max(p0);
end

if numel(rp)~=n,
    rp = randperm(n)';
end

plotCellData(G, rp(p0), 'EdgeColor','k','EdgeAlpha',.05);
outlineCoarseGrid(G, p0)
plotWell(G, W, 'height', 30,'FontSize',12)
axis tight off, view(250, 50),
colormap(.8*tatarizeMap(n)+.2*ones(n,3)); caxis([1 n]);

bv = accumarray(p0,G.cells.volumes); bv=bv/min(bv);
axes('Position',[.02 .82 .35 .16]);
bar(bv,1,'EdgeColor','none');
box on, 
set(gca,'YAxisLocation','right','Fontsize',14,'XLim',[1 max(p0)],...
    'YLim',[0 1.15*max(bv)]);

end

