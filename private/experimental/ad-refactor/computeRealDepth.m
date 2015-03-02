function depth = computeRealDepth(Gt, slope, slopedir, h)
    % Compute depth, taking the inclination of the aquifer into account
    if slope == 0
        % no correction for angle needs to be made
        depth = Gt.cells.z + h;
    else
        ref_cell = [1, 1]; % @@ For now, we use this corner cell as the reference cell
        ref_ix = sub2ind(Gt.cartDims, ref_cell(1), ref_cell(2));
        shift = bsxfun(@minus, Gt.cells.centroids, Gt.cells.centroids(ref_ix, :));
        shift = shift * slopedir(:);

        depth = Gt.cells.z(ref_ix) + ...
            (Gt.cells.z - Gt.cells.z(ref_ix) + h) * cos(slope) - ...
            shift * sin(slope);
    end
end

