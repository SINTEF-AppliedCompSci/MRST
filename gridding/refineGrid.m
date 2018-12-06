function [G, mappings, partition] = refineGrid(G, GC, GF, cells)

    % Mappings (fine grid) <-> (old grid)
    partition = G.partition;
    fine2old  = partition;
    
    GF.cells.num;
    old2fine  = nan(G.cells.num,1);
    old2fine(partition) = (1:GF.cells.num);
    
    fine2coarse = GC.partition;
    old2coarse = fine2coarse(old2fine);
    cells = old2coarse(cells);
    
    partition     = GC.partition;
    ix = any(partition == cells',2);
    
    % Assign new numbers to non-refined cells
    partition(ix) = max(partition) + 1;
    map = nan(max(partition), 1);
    up = unique(partition);
    map(up) = 1:numel(up);
    
    partition     = map(partition);
    partition(ix) = (1:nnz(ix)) + numel(up) - 1;
    
    G = generateCoarseGrid(GF, partition);
    G = coarsenGeometry(G);
%     G = storeInteractionRegionCart(G);
%     
    new2fine = nan(G.cells.num,1);
    new2fine(partition) = (1:GF.cells.num);
    fine2new = partition;
    
    new2old = fine2old(new2fine);
    old2new = fine2new(old2fine);
    
    G.refined = false(G.cells.num,1);
    G.refined(old2new(cells)) = true;
    
    mappings = struct('new2fine', new2fine, ...
                      'fine2new', fine2new, ...
                      'old2fine', old2fine, ...
                      'fine2old', fine2old, ...
                      'new2old' , new2old , ...
                      'old2new' , old2new );

end