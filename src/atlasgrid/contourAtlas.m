function contourAtlas(info, varargin)
% Plot contour lines in 3D for height data
%
% SYNOPSIS:
%       contourAtlas(dataset)
%       contourAtlas(dataset, N)
%       contourAtlas(dataset, N, linewidth)
%
% PARAMETERS:
%   dataset   - Dataset as defined by the second output of getAtlasGrid.
%
%   N         - (OPTIONAL) Number of isolines. Default: 10
%
%   linewidth - (OPTIONAL) Width of lines. Default: 1
%
% RETURNS:
%   Nothing. Will only produce a plot.

%{
#COPYRIGHT#
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
    [C h] = contour3(X,Y,d, N, '-'); %#ok
    set(h,'LineWidth', 1);
    
    data = get(h, 'UserData');
    dpts = unique([data{:}]);

    colors = flipud(jet(N+1));
    for i = 1:numel(h)
        set(h(i), 'LineWidth', l)
        set(h(i), 'Color', colors(dpts == get(h(i), 'UserData'), :));
    end
end
