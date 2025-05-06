function [hc] = plotReservoirModel(G, p, W, prms)
%Make a plot of a reservoir model

if isempty(p)
    plotGrid(G,prms.pargs{:});
else
    plotCellData(G,p,prms.pargs{:});
end
if ~isempty(W)
    if isfield(prms,'wargs')
        plotWell(G, W, prms.wargs{:});
    else
        plotWell(G, W);
    end
end
set(gca,'Clipping', 'off');
drawnow; axis tight off
if isfield(prms,'view'),    view(prms.view); end
if isfield(prms,'dataasp'), set(gca,'DataAspect', prms.dataasp); end
if isfield(prms,'zoom'),    zoom(prms.zoom); end
if isfield(prms,'dolly')
    camdolly(prms.dolly(1), prms.dolly(2), prms.dolly(3));
end
if all(isfield(prms, {'cs', 'unit', 'cb'}))
    caxis(log10([min(prms.cs) max(prms.cs)]*prms.unit));
    hc = colorbarHist(p(~isinf(p)), caxis, prms.cb,100);
    set(hc, 'YTick', 0.5, 'YTickLabel','mD', ...
        'XTick', log10(prms.cs*prms.unit), ...
        'XTickLabel', num2str(prms.cs'));
elseif isfield(prms,'cb')
    hc = colorbar(prms.cb);
end
end