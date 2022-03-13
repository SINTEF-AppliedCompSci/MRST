function varargout = contourAtlas(info, varargin)
% Plot contour lines in 3D for height data
%
% SYNOPSIS:
%       contourAtlas(dataset)
%       contourAtlas(dataset, N)
%       contourAtlas(dataset, N, linewidth)
%   h = contourAtlas(dataset, N, linewidth, color)
%
% PARAMETERS:
%   dataset   - Dataset as defined by the second output of getAtlasGrid.
%
%   N         - (OPTIONAL) Number of isolines. Default: 10
%
%   linewidth - (OPTIONAL) Width of lines. Default: 1
%
%   color     - (OPTIONAL) Set one color for all contour lines.
%               Default: color according to depth
%
% RETURNS:
%   h         - handle to graphics produced by the routine

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
    if nargin == 1
        N = 10;
    else
        N = varargin{1};
    end
    
    if nargin < 3
        l = 1;
    else
        l = varargin{2};
    end
    
    if nargin < 4
       color = [];
    else
       color = varargin{3};
    end
    
    hold on
    set(gca, 'ZDir', 'reverse')

    m = info.meta;
    d = info.data;

    h = m.cellsize;
    xl = m.xllcorner;
    yl = m.yllcorner;
    dims = [m.nrows m.ncols];
    [X, Y]  = ndgrid(linspace(xl, (dims(1) - 1)*h + xl, dims(1) ), ...
                     linspace(yl, (dims(2) - 1)*h + yl, dims(2) ));
    
    % Lineplot to avoid colormap messup
    [C, h] = contour3(X,Y,d, N, '-');                           %#ok<ASGLU>
    set(h,'LineWidth', 1);
    
    if isempty(color)
       data = get(h, 'UserData');

       if ~ isempty(data),
          if iscell(data),
             data = [ data{:} ];
          end

          dpts = unique(data);
          colors = flipud(jet(N+1));
          for i = 1:numel(h)
             set(h(i), 'LineWidth', l)
             set(h(i), 'Color', colors(dpts == get(h(i), 'UserData'), :));
          end
       end
    else
       set(h,'Color',color);
    end

    if nargout > 0,
       varargout{1} = h;
    end
end
