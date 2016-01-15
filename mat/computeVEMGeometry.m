function G = computeVEMGeometry(G)

    edgeVec   = G.nodes.coords(G.edges.nodes(2:2:end),:) -  ...
                 G.nodes.coords(G.edges.nodes(1:2:end-1),:);
    lengths   = sqrt(sum(edgeVec.^2,2));
    centroids = (G.nodes.coords(G.edges.nodes(2:2:end),:) +  ...
                 G.nodes.coords(G.edges.nodes(1:2:end-1),:))./2;
    
    faceNormals = G.faces.normals;
    edgeNum = mcolon(G.faces.edgePos(1:end-1),G.faces.edgePos(2:end)-1);
    edges = G.faces.edges(edgeNum);
    signs = G.faces.edgeSign(edgeNum);
    nodeNum = mcolon(G.edges.nodePos(edges),G.edges.nodePos(edges+1)-1);
    nodes = G.edges.nodes(nodeNum);
    edgeVec = G.nodes.coords(nodes(2:2:end),:)-G.nodes.coords(nodes(1:2:end-1),:);
    edgeVec = edgeVec.*repmat(signs,1,3);
    normals = cross(edgeVec, rldecode(faceNormals, diff(G.faces.edgePos), 1));
    normals = normals./repmat(sqrt(sum(normals.^2,2)),1,3);
    
    cellDiameters = zeros(G.cells.num,1);
    for i = 1:G.cells.num
        nodeNum = G.cells.nodePos(i):G.cells.nodePos(i+1)-1;
        nodes = G.cells.nodes(nodeNum);
        X = G.nodes.coords(nodes,:);
        cellDiameters(i) = cellDiameter(X);
    end
    faceDiameters = zeros(G.faces.num,1);
    for i = 1:G.faces.num
        nodeNum = G.faces.nodePos(i):G.faces.nodePos(i+1)-1;
        nodes = G.faces.nodes(nodeNum);
        X = G.nodes.coords(nodes,:);
        faceDiameters(i) = cellDiameter(X);
    end
    
    G.edges.('lengths')     = lengths;
    G.edges.('centroids')   = centroids;
    G.faces.('edgeNormals') = normals;
    G.cells.('diameters')   = cellDiameters;
    G.faces.('diameters')   = faceDiameters;
    
    [intD, intB] = faceInt(G);
    G.faces.('faceInt') = {intD, intB};
end

