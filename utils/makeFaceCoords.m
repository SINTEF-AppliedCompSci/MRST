function G = makeFaceCoords(G)

    nfe   = diff(G.faces.edgePos);
    edges = G.faces.edges(cumsum(nfe));
    nodes = G.edges.nodes(mcolon(G.edges.nodePos(edges), G.edges.nodePos(edges+1)-1));
    
    v1    = G.nodes.coords(nodes(2:2:end),:) - G.nodes.coords(nodes(1:2:end-1),:);
    v1    = v1./sqrt(sum(v1.^2,2));
    n     = G.faces.normals./G.faces.areas;
    v2    = cross(v1, n, 2);
    
    G.faces.v = {v1, v2};
    
end