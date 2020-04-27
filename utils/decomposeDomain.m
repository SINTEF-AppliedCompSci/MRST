function partition = decomposeDomain(model, W, varargin)
    opt = struct('type'          , 'metis', ...
                 'pdims'         , []     , ...
                 'merge'         , true   , ...
                 'partitionWells', true   , ...
                 'rWell'         , 50     , ...
                 'minBlockSize'  , -inf   );
                 
    opt = merge_options(opt, varargin{:});
    if opt.minBlockSize < 0
        opt.minBlockSize = model.G.cells.num/20;
    end
    % Partition into near-well regions and reservorir
    if opt.partitionWells
        [pw, is_well] = makeWellPartition(model, W, opt);
    else
        pw      = zeros(model.G.cells.num,1);
        is_well = false(model.G.cells.num,1);
    end
    % Partition everything ellse
    [pg, is_grid] = makeGridPartition(model, is_well, opt);
    % Combine the partitions
    partition = pw;
    partition(is_grid) = pg + max(partition);
    partition = compressPartition(partition);
    % Merge too large blocks
    ncmin = opt.minBlockSize;
    if opt.partitionWells
        ncmin = min(accumarray(pw, 1));
    end
    if opt.merge
        partition = mergePartition(model, partition, ncmin);
    end
    nc = accumarray(partition,1);
    partition(nc(partition) < 0.1*ncmin) = max(partition);
    partition = compressPartition(partition);
end

%-----------------------------------------------------------------%
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

%-----------------------------------------------------------------%
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
            pdims = opt.pdims;
            if isempty(pdims)
                n = ceil(sqrt(nb));
                pdims = ones(1, G.griddim);
                pdims(1:2) = n;
            end
            if any(strcmpi(G, 'cartGrid'))
                pg = partitionUI(Gr, pdims);
            else
                pg = partitionCartGrid(pdims*50, pdims);
                pg = sampleFromBox(Gr, reshape(pg, pdims*50));
            end
    end

end

%-----------------------------------------------------------------%
function partition = mergePartition(model, partition, minBlockSize)
    T = getFaceTransmissibility(model.G, model.rock);
    partition = mergeBlocksByConnections(model.G, partition, T, minBlockSize);
end