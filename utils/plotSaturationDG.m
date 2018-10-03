function [h, saturation, coords, keep, n] = plotSaturationDG(disc, state, varargin)

    opt = struct('n'      , 100, ...
                 'phaseNo', 1, ...
                 'plot'   , true, ...
                 'plot1d' , false);
    [opt, extra] = merge_options(opt, varargin{:});
            
    G = disc.G;
    assert(G.griddim == 2);

    n = opt.n;
    xmax = max(G.nodes.coords, [], 1);
    xmin = min(G.nodes.coords, [], 1);
    x = linspace(xmin(1), xmax(1), n)';
    if opt.plot1d
        y = repmat((xmin(2) + xmax(2))/2, size(x,1), 1);
        x = [x,y];
    else
        y = linspace(xmin(2), xmax(2), n);
        [x0,y0] = ndgrid(x,y);
        x = [x0(:), y0(:)];
    end
    nPts = size(x,1);
    
    cellNo = ones(nPts,1);
    keep = false(nPts,1);
    for cNo = 1:G.cells.num
        nodes = getCellNodes(cNo, G);
        xn = G.nodes.coords(nodes,:);
        [in, on] = inpolygon(x(:,1), x(:,2), xn(:,1), xn(:,2));
        cellNo(in | on) = cNo;
        keep = keep | in | on;
    end
    
    xt = disc.transformCoords(x, cellNo);
    
    saturation = @(state) disc.evaluateSaturation(xt, cellNo, state.sdof(:,opt.phaseNo), state);
    if opt.plot1d
        coords = {x};
    else
        coords = {x0, y0};
    end
    
    h = [];
    if opt.plot
        if opt.plot1d
            s = saturation(state);
            h = plot(x(:,1), s, extra{:});
        else
            s = saturation(state);
            s(~keep) = nan;
            s = reshape(s', [n, n]);
            h = surf(x0, y0, s', extra{:});
        end
    end

end

function nodes = getCellNodes(cell, G)
    
    faces = G.cells.faces(G.cells.facePos(cell):G.cells.facePos(cell+1)-1);
    nodes = G.faces.nodes(mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces+1)-1));
    nodes = reshape(nodes, 2, [])';
    swap = G.faces.neighbors(faces,1) ~= cell;
    nodes(swap,:) = nodes(swap, [2,1]); nodes = nodes(:,1);
    
end
    