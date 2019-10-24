function [h, saturation] = plotSaturationDG(disc, state, varargin)

    opt = struct('n'      , 100  , ...
                 'phaseNo', 1    , ...
                 'plot1d' , false, ...
                 'coords' , []   , ...
                 'contour', false, ...
                 'makePlot', true);
    [opt, extra] = merge_options(opt, varargin{:});
    
    coords = opt.coords;
    if isempty(coords)
        coords = getPlotCoordinates(disc.G, varargin{:});
    end
    
    x     = coords.points;
    keep  = coords.keep;
    cells = coords.cells;
    n     = coords.n;
    xs    = cell(1, disc.G.griddim);
    for d = 1:disc.G.griddim
        xs{d} = reshape(x(:,d), n);
    end
    saturation = @(state) disc.evaluateDGVariable(x, cells, state, state.sdof(:,opt.phaseNo));
    h = [];
    if opt.makePlot && ~isempty(state)
        s = saturation(state);
        if opt.plot1d
            h = plot(x(:,1), s, extra{:});
        elseif disc.G.griddim == 2
            s(~keep) = nan;
            s        = reshape(s, n);
            if opt.contour
                h = contourf(xs{1}, xs{2}, s, extra{:});
            else
                h = surf(xs{1}, xs{2}, s, extra{:});
            end
        else
            s(~keep) = nan;
            s        = reshape(s, n);
            h = patch(isosurface(xs{1}, xs{2}, xs{3}, s, 0.6));
            h.FaceColor = 'blue';
            h.EdgeColor = 'none';
            lighting gouraud
            axis equal
            axis([min(x(:,1)), max(x(:,1)), min(x(:,2)), max(x(:,2)), min(x(:,3)), max(x(:,3))]);
        end
    end

end