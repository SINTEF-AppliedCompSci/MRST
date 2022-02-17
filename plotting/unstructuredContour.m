function varargout = unstructuredContour(G, data, varargin)
%Plot data on 2D grid as a surface plot.
%
% SYNOPSIS:
%       unstructuredContour(G, data)
%       unstructuredContour(G, data, n)
%       unstructuredContour(G, data, v)
%       unstructuredContour(G, 'pn1', pv1, ...)
%   h = unstructuredContour(...)
%
% PARAMETERS:
%   G       - Grid data structure.
%
%   data    - Scalar data with which to colour the grid. One scalar,
%             indexed colour value for each cell in the grid or one
%             TrueColor value for each cell. The number of data values
%             can be one of the following:
%               1) G.nodes.num
%               2) G.faces.num
%               3) G.cells.num
%             In cases 2 and 3, the data is linearly interpolated to the
%             node coordinates.
%
%   n/v     - If scalar n, the function draws n coutour lines. If vector,
%             the function draws contour lines for each value in v.
%             OPTIONAL. If empty or not provided, contour/contourf will
%             pick suitable contour values.
%
% OPTIONAL ARGUMENTS:
%   'cartDims' - Countours are plotted on structured grid with cartDims(1)
%                points in x each direction and cartDims(2) points in y
%                direction. If not provided, the function infers suitable
%                values from the grid, but with a maximum of 1000.
%
%   'fill'     - If true, make filled contour plot (using contourf).
%                Default is true.
%
%   'Any'      - Additional keyword arguments will be passed directly on to
%                function `patch` meaning all properties supported by
%                `patch` are valid.
%
% RETURNS:
%   c   - Contour matrix (see conoutrc)
%
%   h   - Handle to resulting contour object. The contour object is added
%         directly to the current `axes` object (`gca`).
%         OPTIONAL. Only returned if specifically requested.
%
% EXAMPLE:
%   % Given a grid 'G' and a reservoir solution structure 'resSol' returned
%   % from, e.g., function 'incompTPFA', plot the cell pressure in bar:
%
%      figure, unstructuredContour(G, convertTo(resSol.pressure, barsa()));
%
% SEE ALSO:
%   `plotCellData`, `patch`

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

    % Optional input arguments
    %---------------------------------------------------------------------%
    v = [];
    if nargin > 2 && ~ischar(varargin{1})
        v = varargin{1}; varargin = varargin(2:end);
    end
    opt = struct('cartDims'     , []      , ...
                 'fill'         , false   , ...
                 'extrapolation', 'linear');
    [opt, extra] = merge_options(opt, varargin{:});
    % Check that we have a 2D grid
    assert(G.griddim == 2, 'unstructuredContour only supports 2D grid');
    %---------------------------------------------------------------------%

    % Make contour plot
    %---------------------------------------------------------------------%
    % Make structured grid and interpolate the unstructured data onto it
    [xq, yq, opt] = makeNDGrid(G, opt);
    data          = interpolateData(G, data, [], [xq(:), yq(:)], ...
                                    'extrapolation', opt.extrapolation);
    % Plot
    if isempty(v), v = 10; end
    if opt.fill
        cfun = @(varargin) contourf(varargin{:});
    else
        cfun = @(varargin) contour(varargin{:});
    end
    [c, h] = cfun(xq, yq, reshape(data, opt.cartDims), v, extra{:});
    
    %---------------------------------------------------------------------%
    
    % Return handle to patch object
    %---------------------------------------------------------------------%
    if nargout > 0, varargout{1} = c; end
    if nargout > 1, varargout{2} = h; end
    %---------------------------------------------------------------------%

end

function [x, y, opt] = makeNDGrid(G, opt)
% Make structured grid on which to plot the contours

    xmin = min(G.nodes.coords);
    xmax = max(G.nodes.coords);
    if isempty(opt.cartDims)
        [dx, dy] = cellDims(G, 1:G.cells.num);
        dx = min([dx, dy],[],1);
        DX = xmax - xmin;
        opt.cartDims = ceil(DX./dx);
    end
    
    nMax = 1000;
    if any(opt.cartDims > nMax)
        warning(['Number of interpolation points exceeds maximum. ', ...
                 'Reducing to 1000 points in each direction. ', ...
                 'Use optional argument cartDims to set this manually']);
        opt.cartDims = min(opt.cartDims, nMax);
    end
    
    [x,y] = ndgrid(linspace(xmin(1), xmax(1), opt.cartDims(1)), ...
                   linspace(xmin(2), xmax(2), opt.cartDims(2)));
    
end