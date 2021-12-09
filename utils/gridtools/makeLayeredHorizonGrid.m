function [G, cellmap] = makeLayeredHorizonGrid(G2D, horizons, varargin)
%Build 3D grid that adapts to a set of horizons from a 2D grid by making a
%layered grid and interpolating z coordinates between horizons
%
% SYNOPSIS:
%   [G, cellmap] = makeLayeredHorizonGrid(G2D, horizons)
%   [G, cellmap] = makeLayeredHorizonGrid(G2D, horizons, 'pn', pv)
%
% PARAMETERS:
%  G2D      - A 2D grid to be extruded to 3D
%  horizons - A cell array of structures describing the surfaces that make
%             up the individual horizons. There must at least be two such
%             horizons. The structure should contain three fields (x,y,z)
%             that give the coordinates of the respective horizons. The x
%             and y coordinates are assumed to lie on a rectilinear grid
%             (i.e., what MATLAB refers to as a meshgrid), but the number
%             of mesh points need not be the same in the different surfaces.
%
%  'layers' - A vector specifying the number of grid layers to be inserted
%             in between each pair of horizon surfaces
%
%  'repairFunction'  - Function handle to repair inter-layer collisions.
%                      Can be e.g. @max or @min to take either the top or
%                      bottom layer for a segment.
%
%  'repairFunction2' - Function handle to repair inter-layer collisions.
%                      Can be e.g. @max or @min to take either the top or
%                      bottom layer for a segment.
%
%  'method' - Method to use for interpolating between points. Can use any
%             of the choices available for griddedInterpolant. Defaults to
%             linear.
%
%  'extrapMethod' - Extrapolation method. Defaults to 'none'.
%
% DESCRIPTION:
%
% RETURNS:
%   G       - A valid MRST grid structure
%   cellmap - A map from cells in the grid G3D with collapsed cells due to
%             collapsing horizons, (with G2D.cells.num*sum(layers) numer of
%             cells), to cells in the returned grid where these are
%             removed. Specifically, `cellmap(i)` is the cell ID of `G3D`
%             that corresponds to cell `i` in `G`.
%
% EXAMPLES:
%   [y,x,z]  = peaks(30);
%   horizons = {struct('x',x,'y',y,'z',z),struct('x',x,'y',y,'z',z+5)};
%   G2D      = tensorGrid(x(:,1), y(1,:)');
%   G        = makeLayeredHorizonGrid(G2D, horizons,'layers', 2);
%   figure, plotGrid(G); view(3); axis tight
%
%   horizons = {struct('x',x,'y',y,'z',z),struct('x',x,'y',2*y+1,'z',z+10)};
%   mrstModule add upr
%   G2D = pebiGrid2D(6/25, [6,6], 'polyBdr', [-3,-3; 3,-3; 3,3; -3,3]);
%   G   = makeLayeredHorizonGrid(G2D, horizons,'layers', 2);
%   figure, plotGrid(G); view(-140,30); axis tight
%
% SEE ALSO:
%   `convertHorizonsToGrid`

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

    opt = struct('layers'         , []      , ... % Number of cells in each layer
                 'method'         , 'linear', ... % Interpolation method
                 'extrapMethod'   , 'none'  , ... % Extrapolation method
                 'doRepair'       , true    , ... % Repair horizon intersections
                 'repairFunction' , []      , ...
                 'repairFunction2', []      ); 
    opt = merge_options(opt, varargin{:});
    
    % Process layer input
    nHorizons = numel(horizons);
    if isempty(opt.layers)
        opt.layers = ones(nHorizons-1, 1);
    elseif numel(opt.layers) == 1
        opt.layers = repmat(opt.layers, nHorizons-1, 1);
    end
    assert(numel(opt.layers) == nHorizons - 1,...
        'Please supply one layer count per layer (= number of horizons-1)');
    % Go through each layer between a pair of horizons and add 
    zHorizons = cellfun(@(h) interpolate(h, G2D, opt), horizons, 'UniformOutput', false);
    if opt.doRepair
        zHorizons = repair(zHorizons, opt.repairFunction, 2:nHorizons, -1);
        zHorizons = repair(zHorizons, opt.repairFunction2, nHorizons-1:-1:1, 1);
    end
    % Make layered grid
    G       = makeLayeredGrid(G2D, sum(opt.layers));
    if isfield(G.faces, 'tag')
        G.faces = rmfield(G.faces, 'tag');
    end
    layer = reshape(rldecode(repmat((1:numel(opt.layers)), G2D.cells.num, 1), opt.layers, 2), [], 1);
    % Set top horizon
    zPrev      = zHorizons{1};
    layerIndex = 1;
    n2d        = G2D.nodes.num;
    nix        = (1:n2d);
    G.nodes.coords(nix,3) = zPrev;
    % Loop through layers and interpolate z coordinates to horizons
    for horizonIndex = 2:nHorizons
        % Find cell layer thickness
        zNext  = zHorizons{horizonIndex};
        nLocal = opt.layers(horizonIndex-1);
        dz     = (zNext - zPrev)/nLocal;
        checkDeltaZ(dz);
        for ix = 1:nLocal
            % Find node index
            nix = (1:n2d) + n2d*layerIndex;
            % Adjust z coordinate
            G.nodes.coords(nix,3) = zPrev + ix*dz;
            % Update layer index
            layerIndex = layerIndex + 1;
        end
        zPrev = zNext;
    end
    % Postprocess to remove collapsed cells
    tol           = (max(zNext) - min(zHorizons{1}))*1e-16;
    [G, cellmap]  = removeShortEdges(G, tol);
    % Preserve G.cells.layers map
    G.cells.layer = layer(cellmap);
end

%-------------------------------------------------------------------------%
function z = repair(z, fn, iter, offset)
    if isempty(fn)
        fn = @(layer, other) layer;
    end
    for i = iter
        self  = z{i};
        other = z{i+offset};
        if offset < 0
            delay       = isinf(self);
            self(delay) = other(delay);
        end
        bad  = isnan(self);
        self = fn(self, other);
        self(bad) = nan;  % Preserve NaN values
        z{i}      = self; % Replace
    end
end

%-------------------------------------------------------------------------%
function F = interpolate(horizon, G, opt)
    X = horizon.x(:);
    Y = horizon.y(:);
    Z = horizon.z(:);
    F = scatteredInterpolant(X, Y, Z, opt.method);
    F = F(G.nodes.coords(:,1), G.nodes.coords(:,2));
end

%-------------------------------------------------------------------------%
function checkDeltaZ(dz)
    dz_finite = dz(:);
    dz_finite = dz_finite(isfinite(dz));
    if ~all(dz_finite >= 0)
        warning('Non-monotone data detected.');
    end
end