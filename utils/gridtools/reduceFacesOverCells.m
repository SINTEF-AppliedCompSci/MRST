function values = reduceFacesOverCells(G, cellNo, neighborship, celldata, fun)
    if nargin < 5
        fun = @sum;
    end

    cellInterface = G.cells.faces(:,1);
    data = celldata(cellNo);
    
    if isfield(G, 'nnc') && isfield(G.nnc, 'cells')
        nnc_cells = [G.nnc.cells(:, 1); G.nnc.cells(:, 2);];
        nnc_faceno = G.faces.num + (1:size(G.nnc.cells, 1)).';
        
        cellInterface = [cellInterface; nnc_faceno; nnc_faceno];
        data = [data; celldata(nnc_cells)]; 
    end
        
    values = accumarray(cellInterface, data, ...
                        [size(neighborship, 1), 1], fun);
    
end
