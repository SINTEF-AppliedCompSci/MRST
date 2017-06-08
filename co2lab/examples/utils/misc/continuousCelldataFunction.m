function [F, mask_fn] = continuousCelldataFunction(Gt, celldata, quickclip)

    % computing mask function
    bnd = Gt.nodes.coords(extractGridBoundaryNodes(Gt),:);
    mask_fn = @(x, y) inpolygon(x, y, bnd(:,1), bnd(:,2));

    % computing interpolating function
    if ~quickclip
        % we will not use boundary notes as NaN values for 'quick' clipping
        bnd = zeros(0,2); 
    end
    
    x = [Gt.cells.centroids(:, 1); bnd(1:end-1,1)];
    y = [Gt.cells.centroids(:, 2); bnd(1:end-1,2)];
    z = [celldata; ones(size(bnd,1)-1,1) * NaN];
    
    F = TriScatteredInterp(x, y, z);
end
