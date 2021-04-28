function [partition, isWell] = makeWellPartition(G, W, varargin)

    opt = struct('rWell'     , 20      , ...
                 'type'      , 'metric', ...
                 'use2d'     , false   , ...
                 'splitWells', true    );
    opt = merge_options(opt, varargin{:});
    % Initialize partition
    p = nan(G.cells.num, numel(W));
    % Get grid connectivity
    N = G.faces.neighbors;
    N = N(all(N>0,2),:);
    M = getConnectivityMatrix(N);
    if opt.use2d, cix = 1:2; else, cix = ':'; end
    x = G.cells.centroids(:, cix);
    for i = 1:numel(W)
        % Add well cells for well i
        cells = false(G.cells.num,1);
        cells(W(i).cells) = true;
        % Make sure we add at least one layer of cells
        cells = cells | M*cells;
        % Add extra cells according to type
        switch opt.type
            case 'metric'
                xw    = x(W(i).cells,:);
                cells = cells | any(pdist2(x, xw) < opt.rWell,2);
            case 'topological'
                for j = 2:opt.rWell
                    cells = cells | M*cells;
                end
        end
        p(:,i) = cells;
    end
    % Make partition
    partition = zeros(G.cells.num, 1);
    if opt.splitWells
        for i = 1:numel(W)
            pi = p(:,i);
            partition(pi > 0) = i;
            for j = 1:numel(W)
                if i == j, continue; end
                pj = p(:,j);
                if any(pi & pj)
                    % Merge intersecting domains
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