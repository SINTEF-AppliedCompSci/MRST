function GC = generateHybridCoarseGrid(G, partition)

    GC = generateCoarseGrid(G, partition);
    GC = coarsenGeometry(GC);
    
    faces = GC.cells.faces(:,1);
    faces = GC.faces.fconn(mcolon(GC.faces.connPos(faces), GC.faces.connPos(faces+1)-1));
    
    nff = diff(GC.faces.connPos);
    ncf = accumarray(rldecode((1:GC.cells.num)', diff(GC.cells.facePos), 1), nff(GC.cells.faces(:,1)));
    GC.cells.faces = faces;
    
    GC.cells.facePos = [0; cumsum(ncf)] + 1;
    
    GC.faces = G.faces;
    GC.faces.neighbors = G.faces.neighbors;
    ix = GC.faces.neighbors ~= 0;
    GC.faces.neighbors(ix) = partition(G.faces.neighbors(ix));
    GC.faces.fconn = (1:GC.faces.num)';
    GC.faces.connPos = (1:GC.faces.num+1)';
    
    GC.edges = G.edges;
    
    GC.nodes = G.nodes;
    
    GC.type = [G.type, { mfilename }];
    
end