function values = reduceFacesOverCells(G, neighborship, cellfacedata, fun)
    if nargin < 5
        fun = @sum;
    end

    cellInterface = G.cells.faces(:,1);
    
    if isfield(G, 'nnc') && isfield(G.nnc, 'cells')
        nnc_cells = [G.nnc.cells(:, 1); G.nnc.cells(:, 2);];
        nnc_faceno = G.faces.num + (1:size(G.nnc.cells, 1)).';
        
        cellInterface = [cellInterface; nnc_faceno; nnc_faceno];
        cellfacedata = [cellfacedata; cellfacedata(nnc_cells)]; 
    end
        
    values = accumarray(cellInterface, cellfacedata, ...
                        [size(neighborship, 1), 1], fun);
    
end
