function unstructuredContour(G, val, varargin)

    fun = scatteredInterpolant(G.cells.centroids(:,1), G.cells.centroids(:,2), val, 'linear', 'none');
    xMax = max(G.nodes.coords);
    xMin = min(G.nodes.coords);
    
    n = 100;
    [x, y] = ndgrid(linspace(xMin(1), xMax(1), n), linspace(xMin(1), xMax(1), n));
    val = fun(x, y);
    
    contour(x, y, val, varargin{:});
    
end
    
