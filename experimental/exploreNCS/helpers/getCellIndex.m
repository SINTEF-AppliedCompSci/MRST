function cellIndex = getCellIndex(Gt, Xcoord, Ycoord)
% Closest cell index of grid Gt corresponding to physical coordinate (X,Y)

    assert(numel(Xcoord) == numel(Ycoord));
    cellIndex = zeros(numel(Xcoord),1);
    
    for i = 1:numel(Xcoord)
        
        dv        = bsxfun(@minus, Gt.cells.centroids(:,1:2), [Xcoord(i), Ycoord(i)]);
        [v, ind]  = min(sum(dv.^2, 2));
        cellIndex(i) = ind; 

    end
    
    % remove any repeated cellIndex @@
    
end
