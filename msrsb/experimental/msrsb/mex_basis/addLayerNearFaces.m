function [p, cells] = addLayerNearFaces(G, p, faces, additive)
    % Add extra partitions near certain fine faces
    if nargin < 4
        additive = true;
    end
    isBF = false(G.faces.num, 1);
    isBF(boundaryFaces(G)) = true;
    
    assert(all(isBF(faces)));
    
    cells = sum(G.faces.neighbors(faces, :), 2);
    
    if additive
        p(cells) = max(p) + p(cells);
    else
        p(cells) = max(p) + 1;
    end
    p = compressPartition(p);
end
