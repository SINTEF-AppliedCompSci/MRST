function mappings = getRefinementMappings(G, GC, GF, cells)
    % G : Current grid
    % GC: Coarse grid
    % GF: Fine grid
    % cells: Cells in G to be refined

    
    if nnz(cells) > 0
        % Assign new numbers to non-refined cells
        if ~islogical(cells)
            c = cells;
            cells = false(G.cells.num,1);
            cells(c) = true;
        end
        cells = GC.partition(cells(G.partition));

        if size(cells,1) > 1
            cells = cells';
        end
        newPartition     = GC.partition;
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
        if isfield(G.cells, 'refined')
            refined      = G.cells.refined;
        else
            refined = [];
        end
    end
    
    new2fine = nan(max(newPartition),1);
    new2fine(newPartition) = (1:GF.cells.num)';
    new2old = G.partition(new2fine);
    
    
    mappings = struct('newPartition', newPartition, ...
                      'oldPartition', G.partition , ...
                      'refined'     , refined     , ...
                      'new2old'     , new2old     );

end