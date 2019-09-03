function [c, h] = unstructuredContour(G, val, nl, varargin)

    opt = struct('useNodeMax', false);
    [opt, extra] = merge_options(opt, varargin{:});

    if nargin == 2
        nl = 10;
    end
    fun = scatteredInterpolant(G.cells.centroids(:,1), G.cells.centroids(:,2), val, 'linear', 'none');
    if opt.useNodeMax
        xMax = max(G.nodes.coords);
        xMin = min(G.nodes.coords); 
    else
        xMax = max(G.cells.centroids);
        xMin = min(G.cells.centroids);
    end
    n = 100;
    n = 20;
    [x, y] = ndgrid(linspace(xMin(1), xMax(1), n), linspace(xMin(2), xMax(2), n));
    val = fun(x, y);
    
    [c, h] = contour(x, y, val, nl, extra{:});
    
end
    
