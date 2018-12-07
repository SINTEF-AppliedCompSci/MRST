function mappings = getRefinementMappings(G, GC, GF, cells)
    % G : Current grid
    % GC: Coarse grid
    % GF: Fine grid
    % cells: Cells in GC to be refined

    
    if ~isempty(cells)
    % Assign new numbers to non-refined cells
        newPartition     = GC.partition;
        if size(cells,1) > 1
            cells = cells';
        end
        ix               = any(newPartition == cells,2);
        newPartition(ix) = max(newPartition) + 1;
        map              = nan(max(newPartition), 1);
        up               = unique(newPartition);
        map(up)          = 1:numel(up);
        newPartition     = map(newPartition);
        newPartition(ix) = (1:nnz(ix)) + numel(up) - 1;

        refined = accumarray(newPartition, ones(GF.cells.num,1)) == 1;
        
    else
        newPartition = G.partition;
        refined      = G.cells.refined;
    end
        
    mappings = struct('newPartition', newPartition   , ...
                      'oldPartition', G.partition, ...
                      'refined'     , refined    );

end