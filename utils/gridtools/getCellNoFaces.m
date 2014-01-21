function [cellNo, cellFaces, isNNC] = getCellNoFaces(G)
    cellNo   = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
    cellFaces = G.cells.faces(:,1);
    
    if isfield(G, 'nnc') && isfield(G.nnc, 'cells')
        nnc_cells = [G.nnc.cells(:, 1); G.nnc.cells(:, 2);];
        nnc_faceno = G.faces.num + (1:size(G.nnc.cells, 1)).';
        cellNo = [cellNo; nnc_cells];
        cellFaces = [cellFaces; nnc_faceno; nnc_faceno];
    end
    isNNC = false(size(cellNo));
    isNNC((end-numel(nnc_cells)+1):end) = true;
end