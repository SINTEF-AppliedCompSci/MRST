function partition = decomposeDomain(model, W, varargin)
    opt = struct('type', 'metis', 'merge', true, 'minBlockSize', -inf, 'rWell', 50);
    opt = merge_options(opt, varargin{:});
    if opt.minBlockSize < 0
        opt.minBlockSize = model.G.cells.num/20;
    end
    % Partition into near-well regions and reservorir
    [pw, is_well] = makeWellPartition(model, W, opt);
    % Partition everything ellse
    [pg, is_grid] = makeGridPartition(model, is_well, opt);
    % Combine the partitions
    partition = pw;
    partition(is_grid) = pg + max(partition);
    partition = compressPartition(partition);
    % Merge too large blocks
    ncw = min(accumarray(pw, 1));
    if opt.merge
        partition = mergePartition(model, partition, ncw);
    end
    nc = accumarray(partition,1);
    partition(nc(partition) < 0.1*ncw) = max(partition);
    partition = compressPartition(partition);
end

function [pw, is_well] = makeWellPartition(model, W, opt)
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

function [pg, is_grid] = makeGridPartition(model, is_well, opt)
    G  = model.G;
    Gr = removeCells(G, is_well);
    is_grid = ~is_well;
    nb = Gr.cells.num/opt.minBlockSize;
    switch opt.type
        case 'metis'
            rock = model.rock;
            rock.perm = rock.perm(is_grid,:);
            T = getFaceTransmissibility(Gr, rock);
            pg = partitionMETIS(Gr, T, nb);
        case 'cart'
            n = ceil(sqrt(nb));
            pg = partitionUI(Gr, [n,n,1]);
    end

end

function partition = mergePartition(model, partition, minBlockSize)
    T = getFaceTransmissibility(model.G, model.rock);
    partition = mergeBlocksByConnections(model.G, partition, T, minBlockSize);
end