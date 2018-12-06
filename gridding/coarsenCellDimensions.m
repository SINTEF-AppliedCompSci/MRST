function G = coarsenCellDimensions(G)

    % Get minimum and maximum cell coordinates
    [~, ix]  = sort(G.partition);
    [~ , nc] = rlencode(G.partition(ix),1);
    [G.cells.xMin, ~           ] = getMinMax(G.parent.cells.xMin(ix,:), nc);
    [~           , G.cells.xMax] = getMinMax(G.parent.cells.xMax(ix,:), nc);
    
    % Calculate bounding box dimensions
    G.cells.dx = G.cells.xMax - G.cells.xMin;
    
    % Get minimum and maximum face coordinates
    nf = diff(G.faces.connPos);
    [G.faces.xMin, ~           ] = getMinMax(G.parent.faces.xMin(G.faces.fconn,:), nf);
    [~           , G.faces.xMax] = getMinMax(G.parent.faces.xMax(G.faces.fconn,:), nf);
    

    G.faces.dx = sqrt(sum((G.faces.xMax - G.faces.xMin).^2, 2));
    
    n = G.faces.normals./sqrt(sum(G.faces.normals.^2,2));
    v = [n(:,2), -n(:,1)];
    x = [G.faces.centroids - G.faces.dx/2.*v;
         G.faces.centroids + G.faces.dx/2.*v];
    ix = [1:G.faces.num; (1:G.faces.num) + G.faces.num];
    x = x(ix(:),:);
     
    G.nodes.coords = x;
    G.faces.nodes = (1:size(x,1))';
    
    stp = 2*(G.griddim-1);
    
    pos = 1:stp:stp*(G.faces.num+1);
    G.faces.nodePos = pos;

    if 0

    clf
    plotGrid(G)
    cNo = 10;
    hold on
    plotGrid(G, cNo);
    plot(G.cells.xMin(cNo, 1), G.cells.xMin(cNo, 2), '.b', 'markerSize', 20);
    plot(G.cells.xMax(cNo, 1), G.cells.xMax(cNo, 2), '.r', 'markerSize', 20);

    f = G.cells.faces(G.cells.facePos(cNo):G.cells.facePos(cNo+1)-1);
    fNo = f(2);
    plot(G.faces.xMin(fNo, 1), G.faces.xMin(fNo, 2), 'o', 'markerSize', 10);
    plot(G.faces.xMax(fNo, 1), G.faces.xMax(fNo, 2), 'o', 'markerSize', 10);
    ff = G.faces.fconn(G.faces.connPos(fNo):G.faces.connPos(fNo+1)-1);
    plotFaces2D(G.parent, ff);
    hold off
    
    clf
    hold on
    plotGrid(G);
    fNo = 40;
    ff = G.faces.fconn(G.faces.connPos(fNo):G.faces.connPos(fNo+1)-1);
    plotFaces2D(G.parent, ff);
    plot(x(fNo,1), x(fNo,2), '.', 'markerSize', 20)
    plot(x(fNo + G.faces.num,1), x(fNo + G.faces.num,2), '.', 'markerSize', 20)
    plot(G.faces.xMax(fNo,1), G.faces.xMax(fNo,2), 'o', 'markerSize', 5)
    plot(G.faces.xMin(fNo,1), G.faces.xMin(fNo,2), 'o', 'markerSize', 5)
    
    hold off
    axis equal tight
    
    end
    
end
