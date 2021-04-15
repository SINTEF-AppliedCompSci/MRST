function grdecl = convertHorizonsToGrid(horizons, varargin)
%Build corner-point grid based on a series of horizons
%
% SYNOPSIS:
%   grdecl = convertHorizonsToGrid(horizons)
%   grdecl = convertHorizonsToGrid(horizons, 'pn', pv)
%
% PARAMETERS:
%  horizons - A cell array of structures describing the surfaces that make
%             up the individual horizons. There must at least be two such
%             horizons. The structure should contain three fields (x,y,z)
%             that give the coordinates of the respective horizons. The x
%             and y coordinates are assumed to lie on a rectilinear grid
%             (i.e., what MATLAB refers to as a meshgrid), but the number
%             of mesh points need not be the same in the different surfaces.
%
%  'dims'   - A 2-vector `[nx, ny]` giving the number of cells in each
%             spatial direction inside the iterpolation area. If this
%             optional parameter is not specified, the number of cells in
%             the interpolation region will be the same as the maximum of
%             mesh points found in each spatial direction for all the input
%             horizons.
%
%  'layers' - A vector specifying the number of grid layers to be inserted
%             in between each pair of horizon surfaces
%
%  'repairFunction' - Function handle to repair inter-layer collisions. Can
%                     be e.g. @max or @min to take either the top or bottom
%                     layer for a segment.
%
%  'method' - Method to use for interpolating between points. Can use any
%             of the choices available for griddedInterpolant. Defaults to
%             linear.
%
%  'extrapMethod' - Extrapolation method. Defaults to 'none'.
%
% DESCRIPTION:
%   The function uses 'interp2' to interpolate between pairs of horizons.
%   The interpolation region is set to be the minimum rectangle that
%   contains the areal bounding boxes of all the horizons, inside which the
%   routine will interpolate on a Cartesian grid of size dims(1) x dims(2).
%   The output grid, however, will only contain cells that are contained
%   inside the maximum areal rectangle that fits inside all the individual
%   areal bounding boxes
%
% RETURNS:
%   grdecl - A GRDECL structure suitable for passing to function
%            `processGRDECL`.
%
% EXAMPLES:
%   [y,x,z]  = peaks(30);
%   horizons = {struct('x',x,'y',y,'z',z),struct('x',x,'y',y,'z',z+5)};
%   grdecl   = convertHorizonsToGrid(horizons,'layers', 2);
%   G        = processGRDECL(grdecl);
%   figure, plotGrid(G); view(3); axis tight
%
%   horizons = {struct('x',x,'y',y,'z',z),struct('x',x,'y',2*y+1,'z',z+10)};
%   G = processGRDECL(convertHorizonsToGrid(horizons,'dims',[20 20]));
%   figure, plotGrid(G); view(-140,30); axis tight
%
% SEE ALSO:
%   `processGRDECL`, `makeModel3`

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    opt = struct('dims',            [],... % X and Y resolution
                 'layers',          [], ... % Number of cells in each layer
                 'method',          'linear', ...
                 'extrapMethod',    'none', ...
                 'doRepair',        true, ...
                 'refinex',         [], ...
                 'refiney',         [], ...
                 'refinez',         [], ...
                 'refineMethod',    'linear', ...
                 'repairFunction',  [],...
                 'repairFunction2', [] ...
                 ); 
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
    [X, Y, Z] = getCoordinates(horizons, ndims, opt);
    ndims = size(X); % Recalculate in case of refinement etc
    dims = ndims - 1;
    grdecl = struct();
    grdecl.cartDims = reshape(dims, 1, []);
    % x, y is uniform for all layers
    x = squeeze(X(:, :, 1));
    y = squeeze(Y(:, :, 2));
    
    % Go through each layer between a pair of horizons and add 
    z_horizons = cellfun(@(h) interpolate(h, x, y, opt), horizons, 'UniformOutput', false);
    if opt.doRepair
        z_horizons = repair(z_horizons, opt.repairFunction, 2:n_horizons, -1);
        z_horizons = repair(z_horizons, opt.repairFunction2, n_horizons-1:-1:1, 1);
    end
    z_prev = z_horizons{1};
    layerIndex = 1;
    
    for horizonIndex = 2:n_horizons
        z_next = z_horizons{horizonIndex};
        nLocal = opt.layers(horizonIndex-1);
        dz = (z_next - z_prev)/nLocal;
        checkDeltaZ(dz);
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

function z = repair(z, fn, iter, offset)
    if isempty(fn)
        fn = @(layer, other) layer;
    end
    for i = iter
        self = z{i};
        other = z{i+offset};
        if offset < 0
            delay = isinf(self);
            self(delay) = other(delay);
        end
        bad = isnan(self);
        self = fn(self, other);
        % Preserve NaN values
        self(bad) = nan;
        z{i} = self; % Replace
    end
end

function F = interpolate(horizon, x, y, opt)
    X = horizon.x;
    Y = horizon.y;
    Z = horizon.z;
    if exist('griddedInterpolant', 'file')
        T = griddedInterpolant(X, Y, Z, opt.method, opt.extrapMethod);
        F = T(x, y);
    else
        F = griddata(X, Y, Z, x, y, opt.method);
    end
end

function checkDeltaZ(dz)
    dz_finite = dz(:);
    dz_finite = dz_finite(isfinite(dz));
    if ~all(dz_finite >= 0)
        warning('Non-monotone data detected.');
    end
end

function [X, Y, Z] = getCoordinates(horizons, ndims, opt)
    x_min_all = cellfun(@(d) min(d.x(:)), horizons);
    x_max_all = cellfun(@(d) max(d.x(:)), horizons);
    
    x_min = min(x_min_all);
    x_max = max(x_max_all);
    
    y_min_all = cellfun(@(d) min(d.y(:)), horizons);
    y_max_all = cellfun(@(d) max(d.y(:)), horizons);
    
    y_min = min(y_min_all);
    y_max = max(y_max_all);
    
    
    % Create grids
    [X, Y, Z]  = ndgrid(discspace1(x_min, x_max, ndims(1), opt.refinex, opt.refineMethod), ...
                        discspace1(y_min, y_max, ndims(2), opt.refiney, opt.refineMethod), ...
                        discspace1(0,  1,        ndims(3), opt.refinez, opt.refineMethod));
end