function monomialVals = monomialValues(G,k)

m3D = retrieve3DMonomials(k);

nodeNum = mcolon(G.cells.nodePos(1:end-1),G.cells.nodePos(2:end)-1);
nodes = G.cells.nodes(nodeNum);
X = G.nodes.coords(nodes,:);
X = bsxfun(@rdivide, X - rldecode(G.cells.centroids, diff(G.cells.nodePos), 1), ...
                         rldecode(G.cells.diameters, diff(G.cells.nodePos), 1));

if k == 2

    edgeNum = mcolon(G.cells.edgePos(1:end-1),G.cells.edgePos(2:end)-1);
    edges = G.cells.edges(edgeNum);
    Ec = G.edges.centroids(edges,:);
    Ec = bsxfun(@rdivide, Ec- rldecode(G.cells.centroids, diff(G.cells.edgePos), 1), ...
                         rldecode(G.cells.diameters, diff(G.cells.edgePos), 1));

end   

if k == 1
    monomialVals = m3D(X);
elseif k == 2
    monomialVals = [m3D(X); m3D(Ec)];
end

end