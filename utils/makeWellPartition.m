function [partition, isWell] = makeWellPartition(G, W, varargin)

    opt = struct('rWell'     , 20      , ...
                 'type'      , 'metric', ...
                 'use2d'     , false   , ...
                 'splitWells', true    );
    opt = merge_options(opt, varargin{:});
    
    p = nan(G.cells.num, numel(W));
    if strcmpi(opt.type, 'topological')
        N = G.faces.neighbors;
        N = N(all(N>0,2),:);
        M = getConnectivityMatrix(N);
    end
    if opt.use2d, cix = 1:2; else, cix = ':'; end
    x = G.cells.centroids(:, cix);
    for i = 1:numel(W)
        switch opt.type
            case 'metric'
                xw    = x(W(i).cells,:);
                cells = any(pdist2(x, xw) < opt.rWell,2);
            case 'topological'
                cells = false(G.cells.num,1);
                cells(W(i).cells) = true;
                for j = 1:opt.rWell
                    cells = cells | M*cells;
                end
        end
        p(:,i) = cells;
    end
    
    partition = zeros(G.cells.num, 1);
    if opt.splitWells
        for i = 1:numel(W)
            pi = p(:,i);
            partition(pi > 0) = i;
            for j = 1:numel(W)
                pj = p(:,j);
                if any(pi & pj)
                    cells = pi | pj;
                    partition(cells) = max(i,j);
                    [p(:,i), p(:,j)] = deal(cells);
                end
            end
        end
        partition = compressPartition(partition)-1;
    else
        partition(any(p,2)) = 1;
    end
    isWell = partition > 0;
    
end