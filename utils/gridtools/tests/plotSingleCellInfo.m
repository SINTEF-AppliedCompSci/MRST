function [] = plotSingleCellInfo(G, c, varargin)
opt = struct('plotFaceInfo', []);
[opt, rest] = merge_options(opt, varargin{:});

fpos = G.cells.facePos([c,(c+1)]);
fix  = G.cells.faces(fpos(1):(fpos(2)-1));

plotInfo = opt.plotFaceInfo;
if ~isempty(plotInfo)
    if ~islogical(plotInfo)
        tmp = false(size(fix));
        tmp(plotInfo) = true;
        plotInfo = tmp;
    end
else
    plotInfo = true(size(fix));
end
    
ax = gca;
if ~ishold(ax)
    hold(ax, 'on');
end

for k = 1:numel(fix)
    plotSingleFaceInfo(G, fix(k), 'fromCell', c, 'plotInfo', plotInfo(k), rest{:})
end
    
tp = get(ax.Children, 'Type');
ix = strcmp(tp, 'text');
uistack(ax.Children(ix), 'top');
end