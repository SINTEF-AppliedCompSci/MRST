function [pw, is_well] = makeWellPartition(model, W, varargin)
    opt = struct('rWell', 20);
    opt = merge_options(opt, varargin{:});
    G = model.G;
    x = G.cells.centroids;
    p = nan(G.cells.num, numel(W));
    for i = 1:numel(W)
        xw = x(W(i).cells,:);
        cells = any(pdist2(x, xw) < opt.rWell,2);
        p(:,i) = cells;
    end
    pw = zeros(G.cells.num, 1);
    for i = 1:numel(W)
        pi = p(:,i);
        for j = 1:numel(W)
            pj = p(:,j);
            if any(pi & pj)
                cells = pi | pj;
                pw(cells) = max(i,j);
                [p(:,i), p(:,j)] = deal(cells);
            end
        end
    end
    is_well = pw > 0;
    pw = compressPartition(pw);
end