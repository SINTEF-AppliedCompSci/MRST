function h = plotDGBasis(G, basis, cellNo, varargin)

    opt          = struct('n', 15);
    [opt, extra] = merge_options(opt, varargin{:});
    
    nodes = G.cells.nodes(G.cells.nodePos(cellNo):G.cells.nodePos(cellNo+1)-1);
    xn    = G.nodes.coords(nodes,:);
    xc = G.cells.centroids(cellNo,:);
    
    xmin = min(xn,[],1);
    xmax = max(xn,[],1);
    
    dx   = max(xmax - xc, xc - xmin)*2;
    
    g = computeGeometry(pebiGrid(max(dx)/opt.n, dx, 'polyBdr', xn));
    x = g.nodes.coords;
    x = (x - G.cells.centroids(cellNo,:))./(dx/2);
    
    z = cell(1, basis.nDof);
    
    xlim = round([xc(1) - dx(1)/2, xc(1) + dx(1)/2, xc(2) - dx(2)/2, xc(2) + dx(2)/2, -1.1, 1.1]*100)/100;
    
    
    faces = boundaryFaces(g);
    for p = 1:basis.nDof
        h(p) = figure();
        plotFaces(g, faces);
        z{p} = basis.psi{p}(x);
        unstructuredCarpetPlot(g, z{p}, 'nodeval', true, extra{:});
        view(3);
        axis tight
        zlim([-1.1, 1.1]);
        pbaspect([1,1,.5]);
        ax = gca;
        ax.XTick = round(linspace(xlim(1), xlim(2), 5)*100)/100;
        ax.YTick = round(linspace(xlim(3), xlim(4), 5)*100)/100;
        caxis([-1,1]);
        colormap('pink')
        box on
        ax.ZDir = 'normal';
    end
    
end