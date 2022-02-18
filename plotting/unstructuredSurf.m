function varargout = unstructuredSurf(G, data, varargin)
%Plot data on 2D grid as a surface plot.
%
% SYNOPSIS:
%       unstructuredSurf(G, data)
%       unstructuredSurf(G, data, 'pn1', pv1, ...)
%       unstructuredSurf(G, data, cells)
%       unstructuredSurf(G, data, cells, 'pn1', pv1, ...)
%   h = unstructuredSurf(...)
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
% OPTIONAL PARAMETERS:
%   'extrapolationMethod' - Method used to extrapolate to nodes on the
%             boundary for case 2 and 3 above. Default is 'linear'
%
%   'Any' -   Additional keyword arguments will be passed directly on to
%             function `patch` meaning all properties supported by `patch`
%             are valid.
%
% RETURNS:
%   h   - Handle to resulting patch object.  The patch object is added
%         directly to the current `axes` object (`gca`).
%         OPTIONAL. Only returned if specifically requested.
%
% EXAMPLE:
%   % Given a grid 'G' and a reservoir solution structure 'resSol' returned
%   % from, e.g., function 'incompTPFA', plot the cell pressure in bar:
%
%      figure, unstructuredSurf(G, convertTo(resSol.pressure, barsa()));
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
    cells = (1:G.cells.num)';
    if nargin > 2 && ~ischar(varargin{1})
        cells    = varargin{1};
        if islogical(cells), cells = find(cells); end
        assert(isnumeric(cells) &&                           ...
               min(cells) >= 1 && max(cells) <= G.cells.num, ...
              'Input argument `cells` must be in the range [1, G.cells.num]');
        varargin = varargin(2:end);
    end
    opt          = struct('extrapolation', 'linear');
    [opt, extra] = merge_options(opt, varargin{:});
    % Check that we have a 2D grid
    assert(G.griddim == 2, 'unstructuredSurf only supports 2D grid');
    %---------------------------------------------------------------------%
    
    % Plot patches
    %---------------------------------------------------------------------%
    % Set node z coordinate equal to data and plot patches
    if numel(data) ~= G.nodes.num
        data = interpolateData(G, data, [], 'nodes', ...
                               'extrapolation', opt.extrapolation);
    end
    G.nodes.coords(:,3) = data;
    h = plotPatches(G, cells, data, extra{:});
    %---------------------------------------------------------------------%
    
    % Return handle to patch object
    %---------------------------------------------------------------------%
    if nargout > 0, varargout{1} = h; end
    %---------------------------------------------------------------------%
    
end