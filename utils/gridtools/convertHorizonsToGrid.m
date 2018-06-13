function grdecl = convertHorizonsToGrid(horizons, varargin)
    opt = struct('dims', [],... % X and Y resolution
                 'layers', []); % Number of cells in each layer
    opt = merge_options(opt, varargin{:});
    n_horizons = numel(horizons);
    
    if isempty(opt.layers)
        opt.layers = ones(n_horizons-1, 1);
    elseif numel(opt.layers) == 1
        opt.layers = repmat(opt.layers, n_horizons-1, 1);
    end
    n_layers = sum(opt.layers);
    assert(numel(opt.layers) == n_horizons - 1,...
        'Please supply one layer count per layer (= number of horizons - 1)');
    
    if isempty(opt.dims)
        nx = cellfun(@(d) size(d.z, 1), horizons);
        ny = cellfun(@(d) size(d.z, 2), horizons);
        ndims = [max(nx), max(ny), n_layers] + 1;
    else
        ndims = [opt.dims(1:2), n_layers] + 1;
    end
    dims = ndims - 1;
    grdecl = struct();
    grdecl.cartDims = reshape(dims, 1, []);
    
    [X, Y, Z] = getCoordinates(horizons, ndims);
    % x, y is uniform for all layers
    x = squeeze(X(:, :, 1));
    y = squeeze(Y(:, :, 2));
    
    % Go through each layer between a pair of horizons and add 
    z_horizons = cellfun(@(h) interpolate(h, x, y), horizons, 'UniformOutput', false);
    z_prev = z_horizons{1};
%     z_prev = interpolate(horizons{1}, x, y);
    layerIndex = 1;
%     Z(:, :, 1) = z_prev;
    
    for horizonIndex = 2:n_horizons
        z_next = z_horizons{horizonIndex};
        nLocal = opt.layers(horizonIndex-1);
        dz = (z_next - z_prev)/nLocal;
        checkDeltaZ(dz)
        
        for ix = 1:nLocal
            Z(:, :, layerIndex) = z_prev + (ix-1)*dz;
            layerIndex = layerIndex + 1;
        end
        z_prev = z_next;
    end
    Z(:, :, end) = z_next;
    Z = sort(Z, 3);

    % Make pillars
    n = prod(ndims(1:2));
    lines = zeros([n, 6]);
    lines(:, [1, 4]) = reshape(X(:,:,[1, end]), [n, 2]);
    lines(:, [2, 5]) = reshape(Y(:,:,[1, end]), [n, 2]);
    lines(:, [3, 6]) = reshape(Z(:,:,[1, end]), [n, 2]);
    grdecl.COORD = reshape(lines.', [], 1);

    % Assign z-coordinates
    % ind(d) == [1, 2, 2, 3, 3, ..., dims(d), dims(d), dims(d)+1]
    ind = @(d) 1 + fix((1 : 2*dims(d)) ./ 2);
    z   = Z(ind(1), ind(2), ind(3));

    grdecl.ZCORN = z(:);

    % Get the cell indices of the Z-corners to find nan values. A cell is
    % only active if all its corners are non-nan (i.e. provided in the
    % original datasets
    cell_index = reshape(1:prod(grdecl.cartDims), grdecl.cartDims);
    cell_index  = rldecode(cell_index, 2, 1);
    cell_index  = rldecode(cell_index, 2, 2);
    cell_index  = rldecode(cell_index, 2, 3);
    
    cell_bad = accumarray(reshape(cell_index, [], 1), grdecl.ZCORN);
    grdecl.ACTNUM = int32(isfinite(cell_bad));
    grdecl.ZCORN(isnan(grdecl.ZCORN)) = inf;
    grdecl.COORD(isnan(grdecl.COORD)) = inf;
end

function F = interpolate(horizon, x, y)
    X = horizon.x;
    Y = horizon.y;
    Z = horizon.z;
    F = interp2(X', Y', Z', x, y);              
end

function checkDeltaZ(dz)
    dz_finite = dz(:);
    dz_finite = dz_finite(isfinite(dz));
    if ~all(dz_finite > 0)
        warning('Non-monotone data detected.');
    end
end

function [X, Y, Z] = getCoordinates(horizons, ndims)
    x_min_all = cellfun(@(d) min(d.x(:)), horizons);
    x_max_all = cellfun(@(d) max(d.x(:)), horizons);
    
    x_min = min(x_min_all);
    x_max = max(x_max_all);
    
    y_min_all = cellfun(@(d) min(d.y(:)), horizons);
    y_max_all = cellfun(@(d) max(d.y(:)), horizons);
    
    y_min = min(y_min_all);
    y_max = max(y_max_all);
    
    
    % Create grids
    [X, Y, Z]  = ndgrid(linspace(x_min, x_max, ndims(1)), ...
                        linspace(y_min, y_max, ndims(2)), ...
                        linspace(0,  1,        ndims(3)));
end