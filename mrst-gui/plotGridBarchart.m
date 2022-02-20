function plotGridBarchart(G, data, varargin)
% Plot barcharts on top of grid for cell data

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
    if mod(numel(varargin), 2) == 1
        cells = varargin{1};
        if islogical(cells)
            cells = find(cells);
        end
        varargin = varargin(2:end);
    else
        cells = (1:G.cells.num) .';
    end

    if isfield(G, 'columns')
        % VE Grid
        cells = G.columns.cells(G.cells.columnPos(cells));
        G = G.parent;
    end

    opt = struct('widthscale',          1, ...
                 'heightscale',         1, ...
                 'edgecolor',           'k', ...
                 'edgealpha',           1, ...
                 'facealpha',           1, ...
                 'sign',                1, ...
                 'scaledata',           true, ...
                 'textformat',          [], ...
                 'fontsize',            12, ...
                 'facecolor',           []);

    opt = merge_options(opt, varargin{:});

    if numel(data) ~= numel(cells)
        assert(numel(data) == G.cells.num);
        data = data(cells);
    end

    % Setup cell -> topface mapping
    faceNum =  rldecode(1 : G.cells.num, diff(G.cells.facePos), 2) .';
    topFaces = zeros(G.cells.num, 1);
    topTag = G.cells.faces(:, 2) == 5;
    topFaces(faceNum(topTag)) = G.cells.faces(topTag, 1);
    topFaces = topFaces(cells);

    % Scale data to unit size
    data = reshape(data, [], 1);
    CData = data;
    if opt.scaledata
        data = data - min(data);
        data = data./max(data);
    end

    % Find the maximum height
    if size(G.cells.centroids, 2) == 3
        % 3D grid
        maxHeight = max(G.faces.centroids(:, 3)) - min(G.faces.centroids(:, 3));
    else
        % Basic 2D grid
        maxHeight = 1;
    end
    maxHeight = maxHeight.*opt.heightscale;
    % scale data into height
    data = opt.sign*data.*maxHeight;

    width = opt.widthscale*sqrt(G.faces.areas(topFaces));
    X = G.cells.centroids(cells, 1);
    Y = G.cells.centroids(cells, 2);
    bottom = G.faces.centroids(topFaces, 3);

    x = fixData(X, X, X, X);
    y = fixData(Y-width/2, Y-width/2, Y+width/2, Y+width/2);
    z = fixData(bottom, bottom - data, bottom - data, bottom);
    d = fixData(CData, CData, CData, CData);

    wasHeld = ishold();
    hold on;
    if isempty(opt.facecolor)
        patch(x, y, z, d,       ...
                                'edgecolor', opt.edgecolor, ...
                                'facealpha', opt.facealpha',...
                                'edgealpha', opt.edgealpha);
    else
        patch(x, y, z, 0*d, 'FaceColor', get_rgb(opt.facecolor),...
                                'edgecolor', opt.edgecolor, ...
                                'facealpha', opt.facealpha',...
                                'edgealpha', opt.edgealpha);
    end
    if ~wasHeld
        hold off;
    end
    set(gca, 'ZDir', 'reverse');

    if ~isempty(opt.textformat)
        for i = 1:numel(X)
            if ischar(opt.textformat)
                name = sprintf(opt.textformat, CData(i));
            else
                name = opt.textformat(CData(i));
            end
            text(X(i), Y(i), bottom(i) - data(i), name, ...
                                        'VerticalAlignment', 'bottom',...
                                        'FontSize', opt.fontsize, ...
                                        'Color', 'r');
        end
    end
end

function x = fixData(a, b, c, d)
    x = [a, b, c, d] .';

end


function rgb = get_rgb(colour)
   if isnumeric(colour)
       rgb = colour;
      return
   end
   switch lower(colour),
      case {'y', 'yellow' }, rgb = [1, 1, 0];
      case {'m', 'magenta'}, rgb = [1, 0, 1];
      case {'c', 'cyan'   }, rgb = [0, 1, 1];
      case {'r', 'red'    }, rgb = [1, 0, 0];
      case {'g', 'green'  }, rgb = [0, 1, 0];
      case {'b', 'blue'   }, rgb = [0, 0, 1];
      case {'w', 'white'  }, rgb = [1, 1, 1];
      case {'k', 'black'  }, rgb = [0, 0, 0];
      otherwise            , rgb = [0, 0, 1]; % Unknown colour -> 'blue'.
   end
end
