function amax=myPlotAllocation(WP, WPc, names, amax)
colormap(.6*jet+.4);
cwp = cumsum(WPc.inj(2).alloc,1); amax = max([sum(cwp,2); amax]);
barh(WPc.inj(2).z, cwp, 'stacked', 'BarWidth', .98, 'EdgeColor','none');
hold on
cwp = cumsum(WP.inj(2).alloc,1);  amax = max([sum(cwp,2); amax]);
if size(cwp,1)<20
    barh(WP.inj(2).z, cwp, 'stacked','BarWidth', .98, 'FaceColor','none');
else
    c = cumsum(cwp,2);
    z =  WP.inj(2).z; dz = min(diff(z));
    n = size(c,2);
    stairs([zeros(1,n); c], repmat(z([1:end end]),1,n), '-k','LineWidth',1);
end
hold off, axis tight
lh=legend(names,'Location','SouthEast');
set(lh,'units','pixels','FontSize',8);
end

