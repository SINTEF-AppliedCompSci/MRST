function [h, saturation] = plotSaturationDG(disc, state, varargin)
    % Visualize dG saturation with higher-order dG basis functions
    opt = struct('n'      , 100  , ...
                 'phaseNo', 1    , ...
                 'plot1d' , false, ...
                 'coords' , []   , ...
                 'contour', false, ...
                 'makePlot', true, ...
                 'pbaspect', [1,1,1]);
    [opt, extra] = merge_options(opt, varargin{:});
    coords = opt.coords;
    if isempty(coords)
        coords = getPlotCoordinates2(disc.G, varargin{:});
    end
    x     = coords.points;
    cells = coords.cells;
    xs    = cell(1, disc.G.griddim);
    saturation = @(state) disc.evaluateDGVariable(x, cells, state, state.sdof(:,opt.phaseNo));
    h = [];
    if opt.makePlot && ~isempty(state)
        s = saturation(state);
        if opt.plot1d
            h = plot(x(:,1), s, extra{:});
        elseif disc.G.griddim == 2
            if opt.contour
                h = contourf(xs{1}, xs{2}, s, extra{:});
            else
                x = [coords.points, s];
                h = patch('faces', coords.faces, 'vertices', x, 'FaceVertexCData', s, 'FaceColor','interp', extra{:});
                zlim([0,1]);
                pbaspect(opt.pbaspect);
            end
        else
            error('Only supported for 1D and 2D')
        end
    end

end