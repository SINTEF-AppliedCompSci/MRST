function showWellCommunication(d, ax, val)
val = val([1:end end],:);
val = double(val(:,[1:end end]));
% scale according to injectors
val = bsxfun(@rdivide, val, sum(val, 2));
val(val==0) = nan;

pcolor(ax, val);
inames= arrayfun(@(x) x.label.String, d.WellPlot.injectors,'UniformOutput',false);
pnames= arrayfun(@(x) x.label.String, d.WellPlot.producers,'UniformOutput',false);
set(ax, 'XTick', 1.5:numel(pnames)+.5, ...
    'XTickLabel', pnames, 'YTick',1.5:numel(inames)+.5, ...
    'YTickLabel', inames, 'XTickLabelRotation',45, ...
    'XAxisLocation','top','Fontsize',8,'YDir','reverse');
axis(ax,'on','tight');
maxVal = max(max(val));
if maxVal == 0 || ~isfinite(maxVal)
    maxVal = 1;
end
caxis(ax, [0 maxVal]);
colorbar(ax,'southoutside')
end