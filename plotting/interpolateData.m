function data = interpolateData(G, data, varargin)
%Interpolate data to corrdinates in a grid
%
% SYNOPSIS:
%   data = interpolateData(G, data)
%   data = interpolateData(G, data, x)
%   data = interpolateData(G, data, x, xq)
%   data = interpolateData(..., 'pn1', vn1, ...)
%
% PARAMETERS:
%   G    - A grid_structure with computed geometry
%
%   data - Dataset to interpolate. Number of rows must match the number
%          of data coordinates x
%
%   x    - Data coordinates. Must equal the number of rows of data
%          OPTIONAL. If empty or not provided, the data coordinates are
%          assumed to be either G.nodes.coordinates, G.faces.centroids or
%          G.cells.centroids, depending on size(data,1)
%
%   xq   - Query coordinates. Can be a set of coordinates, or a string. If
%          string, possible values and corresponding coordinates are
%               * 'nodes' -> G.nodes.coords
%               * 'faces' -> G.faces.centroids
%               * 'cells' -> G.cells.centroids
%          OPTIONAL. If not provided, the default value is 'nodes'
%          
% OPTIONAL PARAMETERS:
%   'interpolation' - Method used to interpolate data. Must be supported by
%                     `scatteredInterpolant`.
%
%   'extrapolation' - Method used to extrapolate data. Must be supported by
%                     `scatteredInterpolant`.
%
% RETURNS:
%   data - Data interpolated to the query points
%
% NOTES:
%   This function is primarily implemented to facilitate visualization, and
%   should be used with care for actual computations.
%
% SEE ALSO:
%   `unstructuredSurf`, unstructuredContour`

    % Process input
    %---------------------------------------------------------------------%
    % Check if data coordinates are provided
    [x, xq] = deal([]);
    if nargin > 2 && ~ischar(varargin{1})
        x = varargin{1}; varargin = varargin(2:end);
    end
    if nargin > 3 && ...
            any(strcmpi(varargin{1}, {'nodes', 'faces', 'cells'})) || ...
            ~ischar(varargin{1})
        xq = varargin{1}; varargin = varargin(2:end);
    end
    opt = struct('interpolation', 'linear', ...
                 'extrapolation', 'linear');
    opt = merge_options(opt, varargin{:});
    % Check that grid has computed geometry
    assert(any(strcmpi(G.type, 'computeGeometry') | ...
               strcmpi(G.type, 'mcomputeGeometry')), ...
               'G must have computed geometry (G = computeGeometry(G))')
    %---------------------------------------------------------------------%
    
    % Interpolate data to query coordinates
    %---------------------------------------------------------------------%
    x    = getDataCoordinates(G, data, x);
    xq   = getQueryCoordinates(G, xq);
    data = interpolateToQuery(x, data, xq, opt);

end

%-------------------------------------------------------------------------%
function x = getDataCoordinates(G, data, x)
% Get data coordinates by looking at the numer of data points

    n = size(data,1);
    if isempty(x)
        switch n
            case G.nodes.num
                x = G.nodes.coords;
            case G.faces.num
                x = G.faces.centroids;
            case G.cells.num
                x = G.cells.centroids;
            otherwise
                error(['Data coordinates x not provided, and size(data,1)', ...
                       ' is not equal to number of nodes, faces or cells ', ...
                       'in the grid. I am not sure how to interpert this']);
        end
    end
    assert(n == size(x,1), ...
        'Number of data points must match number of data coordinates');
    assert(size(x,2) == G.griddim, ...
        'Data coordinates must have G.griddim dimensions')
    
end

%-------------------------------------------------------------------------%
function x = getQueryCoordinates(G, x)
% Get query coordinates from string

    % Default is node coordinates
    if isempty(x), x = 'nodes'; end
    if ischar(x)
        switch x
            case 'nodes'
                x = G.nodes.coords;
            case 'faces'
                x = G.faces.centroids;
            case 'cells'
                x = G.cells.centroids;
            otherwise
                error(['query must be either '              , ...
                       '''nodes'', ''faces'', or ''cells''']);
        end
    end
    assert(size(x,2) == G.griddim, ...
            'Query coordinates must have G.griddim dimensions');
    
end

%-------------------------------------------------------------------------%
function data = interpolateToQuery(x, data, xq, opt)
% Interpolate all columns of data to query

    nc    = size(data,2);
    data0 = data;
    data  = nan(size(xq,1), nc);
    for i = nc
        F = scatteredInterpolant(x, data0(:,i), ...
                                    opt.interpolation, opt.extrapolation);
        data(:,i) = F(xq);
    end
    
end