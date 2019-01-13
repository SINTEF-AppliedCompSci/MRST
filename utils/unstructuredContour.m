function unstructuredContour(G, val, nl, varargin)

    if nargin == 2
        nl = 10;
    end
    fun = scatteredInterpolant(G.cells.centroids(:,1), G.cells.centroids(:,2), val, 'linear', 'none');
    xMax = max(G.cells.centroids);
    xMin = min(G.cells.centroids);
%     xMax = max(G.nodes.coords);
%     xMin = min(G.nodes.coords);
%     
    n = 100;
    [x, y] = ndgrid(linspace(xMin(1), xMax(1), n), linspace(xMin(2), xMax(2), n));
    val = fun(x, y);
    
    contour(x, y, val, nl, varargin{:});
    
end
    
